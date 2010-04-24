(in-package :lla)

(defmacro copy-columns-loop% (&optional (source-element-form 
                                        '(row-major-aref source source-index)))
  "Double loop for copy copying elements columnwise.  Used in
copy-columns etc, see that function for variable names."
  `(iter
     (for (the fixnum col) :from 0 :below ncol)
     (declare (iterate:declare-variables))
     (iter
       (for row :from 0 :below nrow)
       (for source-index
            :from (+ source-offset (cm-index2 source-ld 0 col)))
       (for destination-index
            :from (+ destination-offset (cm-index2 destination-ld 0 col)))
       (declare (iterate:declare-variables)
                (fixnum row source-index destination-index))
       (setf (row-major-aref destination destination-index)
            ,source-element-form))))

(defun copy-columns% (nrow ncol source source-offset source-ld source-type
                      destination destination-offset destination-ld)
  "Fast version of copy-columns for coinciding types."
  (declare (optimize (speed 3) (safety 0))
           (fixnum nrow ncol source-offset source-ld
                   destination-offset destination-ld))
  (expand-for-lla-types (lla-type :prologue (ecase source-type))
    `(,lla-type
      (locally
          (declare (type ,(nv-array-type lla-type) source destination))
        (copy-columns-loop%)))))

(defun copy-columns (nrow ncol source source-offset source-ld source-type
                     destination destination-offset destination-ld
                     &optional (destination-type source-type))
  "Copy elements from one matrix to another, respecting leading
dimensions and offsets.  If types are different, elements are coerced.
Source-type is allowed to be T, in which case elements are always
coerced."
  (if (eq source-type destination-type)
      (copy-columns% nrow ncol
                     source source-offset source-ld source-type
                     destination destination-offset destination-ld)
      (let ((destination-lisp-type (lla-type->lisp-type destination-type)))
        (copy-columns-loop% (coerce (row-major-aref source source-index)
                                    destination-lisp-type))))
  (values))

(defgeneric submatrix (matrix row-start row-end col-start col-end)
  (:documentation "Return submatrix (a DENSE-MATRIX) specified by
row/column start and end indexes (the latter exclusive, with NIL
referring to the last one).  Negative indexes count from the last
column, eg (submatrix matrix 0 0 -1 -1) will drop the last row &
column."))

(defun submatrix-index-calculations% (nrow ncol row-start row-end col-start col-end
                                      &optional (leading-dimension nrow))
  "Return (VALUES NEW-NROW NEW-NCOL OFFSET).  Negative and NIL indexes
are converted according to the SUBMATRIX syntax.  Resulting indexes
are checked for validity."
  (flet ((calc-index (index total start?)
           ;; calculate new index, converting negative specifications
           (cond
             ((null index)
              (if start?
                  (error "~A is not a valid start index" index)
                  total))
             ((minusp index)
              (aprog1 (+ total index)
                (assert (<= 0 it) () "~A~A gives negative index" total index)))
             (t 
              (assert (<= index total) () "index ~A above ~A" index total)
              index))))
    (let* ((row-start (calc-index row-start nrow t))
           (new-nrow (- (calc-index row-end nrow nil) row-start))
           (col-start (calc-index col-start ncol t))
           (new-ncol (- (calc-index col-end ncol nil) col-start))
           (offset (+ (* leading-dimension col-start) row-start)))
      (assert (plusp new-nrow) () "Resulting number of rows is not positive")
      (assert (plusp new-ncol) () "Resulting number of columns is not positive")
      (values new-nrow new-ncol offset))))

(defmethod submatrix ((matrix dense-matrix-like) row-start row-end col-start col-end)
  (declare (optimize (debug 3) (speed 0)))
  (bind (((:accessors-r/o elements lla-type nrow ncol) matrix)
         ((:values new-nrow new-ncol offset)
          (submatrix-index-calculations% nrow ncol row-start row-end
                                         col-start col-end)))
    (set-restricted matrix)
    (let ((result (make-matrix lla-type new-nrow new-ncol)))
      (copy-columns new-nrow new-ncol
                    elements offset nrow lla-type
                    (elements result) 0 new-nrow)
      result)))

(defun subcolumn (matrix col)
  "Return a column from a matrix"
  (copy-nv (submatrix matrix 0 nil col (1+ col)) :copy-p nil))

(defun subrow (matrix row)
  "Return a row from a matrix"
  (copy-nv (submatrix matrix row (1+ row) 0 nil) :copy-p nil))

;;;; transpose
;;;
;;; Methods should take care of returning the correct result type (eg
;;; lower-triangular into upper-triangular, etc).  Helper function
;;; transpose% can be used for implementation.

(defun transpose% (matrix transposed-matrix-kind conjugate-p)
  "Return the transpose of MATRIX, which will be of kind
TRANSPOSED-MATRIX-KIND.  SET-RESTRICTED is *NOT* enforced (so the
caller has to decide whether to enforce it).  Meant to be used as a
helper function, *NOT EXPORTED*."
  (declare (optimize (speed 3) (safety 0)))
  (expand-for-lla-types (lla-type :prologue (ecase (lla-type matrix)))
    (let ((complex-p (lla-complex-p lla-type)))
      `(,lla-type (bind (((:slots-read-only nrow ncol elements) matrix)
                         (transposed (make-nv-elements ,lla-type (length elements))))
                    (declare (fixnum nrow ncol)
                             (type ,(nv-array-type lla-type) elements transposed))
                    (dotimes (col ncol)
                      (declare (fixnum col))
                      (iter
                        (declare (iterate:declare-variables))
                        (for (the fixnum elements-i) :from (cm-index2 nrow 0 col)
                             :below (cm-index2 nrow nrow col))
                        (for (the fixnum transposed-i) :from (cm-index2 ncol col 0) :by ncol)
                        (setf (aref transposed transposed-i)
                              ,(if complex-p
                                   '(let ((elt (aref elements elements-i)))
                                     (if conjugate-p (conjugate elt) elt))
                                   '(aref elements elements-i)))))
                    (make-matrix* ,lla-type ncol nrow transposed 
                                  :kind transposed-matrix-kind))))))

(defgeneric transpose (matrix &optional conjugate-p)
  (:documentation "Return the transpose of a matrix.")
  (:method ((matrix dense-matrix) &optional (conjugate-p t))
    (transpose% matrix :dense conjugate-p))
  (:method ((matrix upper-triangular-matrix) &optional (conjugate-p t))
    (transpose% matrix :lower-triangular conjugate-p))
  (:method ((matrix lower-triangular-matrix) &optional (conjugate-p t))
    (transpose% matrix :upper-triangular conjugate-p))
  (:method ((matrix hermitian-matrix) &optional (conjugate-p t))
    (set-restricted matrix)
    (case (lla-type matrix)
      ((:complex-single :complex-double)
         (transpose% matrix :hermitian conjugate-p))
      (otherwise
         (copy-matrix% matrix :copy-p t)))))

;;;; extraction methods for factorization components and related
;;;; utility functions

(defun matrix-from-first-rows (lla-type vector nrow ncol leading-dimension &optional (kind :dense))
  "Extract & return (as KIND, default is :DENSE) the first NROW rows
of an LEADING-DIMENSIONxNCOL matrix, given as a Lisp vector with
column-major indexing.  NOTE: needed to interface to LAPACK routines
like xGELS."
  (let ((result (make-matrix lla-type nrow ncol :kind kind)))
    (copy-columns% nrow ncol 
                   vector 0 leading-dimension lla-type
                   (elements result) 0 nrow)
    result))


     
;;;; stacking

(defgeneric fill-matrix-elements-using% (object matrix-elements offset
                                                leading-dimension lla-type
                                                horizontal?)
  (:method ((nv numeric-vector) matrix-elements offset leading-dimension
            lla-type horizontal?)
    (bind (((:slots-r/o elements) nv))
      (if horizontal?
          (copy-elements (length elements)
                         elements 0 (lla-type nv)
                         matrix-elements offset lla-type)
          (copy-columns 1 (length elements)
                        elements 0 1 (lla-type nv)
                        matrix-elements offset leading-dimension lla-type))))
  (:method ((diagonal diagonal) matrix-elements offset leading-dimension
            lla-type horizontal?)
    (iter
      (for element :in-vector (elements diagonal) :with-index index)
      (setf (aref matrix-elements (+ offset (cm-index2 leading-dimension index index)))
            (coerce* element lla-type))))
  (:method ((matrix dense-matrix-like) matrix-elements offset leading-dimension
            lla-type horizontal?)
    (bind (((:slots-r/o nrow ncol) matrix))
      (copy-columns nrow ncol
                    (elements matrix) 0 nrow (lla-type matrix)
                    matrix-elements offset leading-dimension)) lla-type))

(defun stack% (objects horizontal?)
  (bind (((:values nrows ncols)
          (iter
            (for object :in objects)
            (for (values nrow ncol) :=
                 (etypecase object
                   (dense-matrix-like
                      (values (nrow object) (ncol object)))
                   (diagonal
                      (let ((l (length (elements object))))
                        (values l l)))
                   (numeric-vector
                      (let ((l (length (elements object))))
                        (if horizontal?
                            (values l 1)
                            (values 1 l))))))
            (collecting nrow :into nrows)
            (collecting ncol :into ncols)
            (finally
             (return (values nrows ncols))))))
    (assert (apply #'= (if horizontal? nrows ncols)) () "Dimensions don't match.")
    (bind ((offset 0)
           (lla-type (apply #'common-target-type (mapcar #'lla-type objects))))
      (if horizontal?
          (bind ((nrow (first nrows))
                 ((:lla-matrix result :elements elements) (make-matrix lla-type nrow (reduce #'+ ncols))))
            (iter
              (for ncol% :in ncols)
              (for object :in objects)
              (fill-matrix-elements-using% object elements offset nrow lla-type t)
              (incf offset (* ncol% nrow)))
            result)
          (bind ((nrow (reduce #'+ nrows))
                 (ncol (first ncols))
                 ((:lla-matrix result :elements elements) (make-matrix lla-type nrow ncol)))
            (iter
              (for nrow% :in nrows)
              (for object :in objects)
              (fill-matrix-elements-using% object elements offset nrow lla-type nil)
              (incf offset nrow%))
            result)))))

(defun stack-horizontally (&rest objects)
  "Stack arguments horizontally, converting to a common type.  A
  vector is interpreted as a column matrix."
  (stack% objects t))

(defun stack-vertically (&rest objects)
  "Stack arguments horizontally, converting to a common type.  A
  vector is interpreted as a column matrix."
  (stack% objects nil))


;;;; identity

(defun eye (n &key (kind :dense) (initial-element 1) (lla-type :double))
  "Return an identity matrix of the given KIND (can also be :diagonal)
and LLA-TYPE.  INITIAL-ELEMENT will be used for the diagonal."
  (if (eq kind :diagonal)
      (make-diagonal lla-type n initial-element)
      (bind (((:lla-matrix eye) (make-matrix lla-type n n :kind kind))
             (initial-element (coerce* initial-element lla-type)))
        (dotimes (i n)
          (setf (eye (eye-index i i)) initial-element))
        eye)))
