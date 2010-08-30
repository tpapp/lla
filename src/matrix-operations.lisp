(in-package :lla)

;;;; extraction methods for factorization components and related
;;;; utility functions

(defun matrix-from-first-rows (vector nrow ncol leading-dimension
                               &optional (kind :dense))
  "Extract & return (as KIND, default is :DENSE) the first NROW rows
of an LEADING-DIMENSIONxNCOL matrix, given as a Lisp vector with
column-major indexing.  NOTE: needed to interface to LAPACK routines
like xGELS."
  (let ((result (make-similar-vector vector (* nrow ncol))))
    (copy-columns nrow ncol 
                  vector 0 leading-dimension
                  result 0 nrow)
    (make-matrix% nrow ncol result :kind kind)))

;;;; stacking

(defgeneric stack-dimensions (object horizontal?)
  (:documentation "Return dimensions when used as a matrix for stacking.
HORIZONTAL? determines the direction for stacking.  Return (values NROW NCOL)")
  (:method ((object dense-matrix-like) horizontal?) 
    (values (nrow object) (ncol object)))
  (:method ((object diagonal) horizontal?)
    (let ((l (length (elements object))))
      (values l l)))
  (:method ((object list) horizontal?)
    (let ((l (length object)))
      (if horizontal?
          (values l 1)
          (values 1 l))))
  (:method ((object vector) horizontal?)
    (let ((l (length object)))
      (if horizontal?
          (values l 1)
          (values 1 l))))
  (:method ((object array) horizontal?)
    (bind (((nrow ncol) (array-dimensions object)))
      (values nrow ncol))))

(defun fill-elements-using-diagonal (elements diagonal-elements start leading-dimension)
  "Copy elements from DIAGONAL to ELEMENTS, using row-major-aref, starting at
START."
  (iter
    (for element :in-vector diagonal-elements)
    (for index :from start :by (1+ leading-dimension))
    (setf (row-major-aref elements index) element)))

(defgeneric fill-stacked-elements (target source horizontal? offset)
  (:documentation "Fill elements of TARGET using SOURCE.  If HORIZONTAL?, OFFSET
gives the starting column, otherwise the starting row.")
  (:method (target (source vector) horizontal? offset)
    (if horizontal?
        (setf (sub target t offset) source)
        (setf (sub target offset t) source)))
  (:method (target (source list) horizontal? offset)
    (if horizontal?
        (setf (sub target t offset) source)
        (setf (sub target offset t) source)))
  (:method (target (source array) horizontal? offset)
    (bind (((nrow ncol) (array-dimensions source)))
      (if horizontal?
          (setf (sub target t (si offset (+ offset ncol))) source)
          (setf (sub target (si offset (+ offset nrow)) t) source))))
  (:method (target (source dense-matrix-like) horizontal? offset)
    (bind (((:slots-r/o nrow ncol) source))
      (if horizontal?
          (setf (sub target t (si offset (+ offset ncol))) source)
          (setf (sub target (si offset (+ offset nrow)) t) source))))
  (:method (target (source diagonal) horizontal? offset)
    ;; !! could be made more efficient by using FILL-ELEMENTS-USING-DIAGONAL,
    ;; !! that would also require that we initialize with 0's
    (fill-stacked-elements target (as-array source) horizontal? offset)))

;; (defgeneric fill-matrix-elements-using% (object matrix-elements offset
;;                                                 leading-dimension
;;                                                 horizontal?)
;;   (:documentation "Fill MATRIX-ELEMENTS (starting at OFFSET, with
;;   LEADING-DIMENSION) using OBJECT.  HORIZONTAL? determines the
;;   direction for stacking.  Return no values.")
;;   (:method ((vector vector) matrix-elements offset leading-dimension
;;             horizontal?)
;;     (if horizontal?
;;         (copy-elements vector 0 
;;                        matrix-elements offset (length vector))
;;         (copy-columns 1 (length vector)
;;                       vector 0 1
;;                       matrix-elements offset leading-dimension)))
;;   (:method ((diagonal diagonal) matrix-elements offset
;;             leading-dimension horizontal?)
;;     (bind (((:slots-r/o elements) diagonal))
;;       ;; dirty hack: we use 1+leading-dimension to "slide" the diagonal
;;       (copy-columns 1 (length elements)
;;                     elements 0 1
;;                     matrix-elements offset (1+ leading-dimension))))
;;   (:method ((matrix dense-matrix-like) matrix-elements offset
;;             leading-dimension horizontal?)
;;     (set-restricted matrix)
;;     (bind (((:slots-r/o elements nrow ncol) matrix))
;;       (copy-columns nrow ncol
;;                     elements 0 nrow
;;                     matrix-elements offset leading-dimension)))
;;   (:method ((matrix array) matrix-elements offset leading-dimension
;;             horizontal?)
;;     (fill-matrix-elements-using% (as-matrix matrix) matrix-elements
;;                                  offset
;;                                  leading-dimension horizontal?)))

(defun stack (target-type direction &rest objects)
  "Stack objects in the given direction (:HORIZONTAL/:H or :VERTICAL/:V).
Return a matrix with the given target type (:DENSE-MATRIX/:MATRIX gives a dense
matrix, :ARRAY a Common Lisp array).  Vectors are automatically treated as
column or row vectors, depending on DIRECTION."
  (declare (optimize debug))
  (bind ((dense-matrix? (ecase target-type
                          ((:matrix :dense-matrix) t)
                          (:array nil)))
         (horizontal? (ecase direction
                        ((:horizontal :h) t)
                        ((:vertical :v) nil)))
         ((:values nrows ncols)
          (iter
            (for object :in objects)
            (for (values nrow ncol) := 
                 (stack-dimensions object horizontal?))
            (collecting nrow :into nrows)
            (collecting ncol :into ncols)
            (finally
             (return (values nrows ncols))))))
    (assert (apply #'= (if horizontal? nrows ncols)) ()
            "Dimensions don't match.")
    (bind ((offset 0)
           (type (reduce #'common-supertype objects
                         :key (lambda (obj)
                                (etypecase obj
                                  (array (array-element-type obj))
                                  (list t)
                                  (elements% (array-element-type (elements obj)))))))
           ((:flet make-result (nrow ncol))
            (if dense-matrix?
                (make-matrix nrow ncol (representable-lla-type type))
                (make-array (list nrow ncol) :element-type type))))
      (if horizontal?
          (bind ((total-ncol (reduce #'+ ncols))
                 (result (make-result (first nrows) total-ncol)))
            (iter
              (for ncol :in ncols)
              (for object :in objects)
              (fill-stacked-elements result object horizontal? offset)
              (incf offset ncol))
            result)
          (bind ((total-nrow (reduce #'+ nrows))
                 (result (make-result total-nrow (first ncols))))
            (iter
              (for nrow :in nrows)
              (for object :in objects)
              (fill-stacked-elements result object horizontal? offset)
              (incf offset nrow))
            result)))))

;;; identity

(defun eye (n &key (kind :dense) (initial-element 1)
            lla-type)
  "Return an identity matrix of the given KIND (can also be :diagonal)
and LLA-TYPE.  INITIAL-ELEMENT will be used for the diagonal."
  (if (eq kind :diagonal)
      (make-diagonal n lla-type initial-element)
      (bind (((:lla-matrix eye) (make-matrix n n lla-type :kind kind))
             (initial-element (coerce* initial-element lla-type)))
        (dotimes (i n)
          (setf (eye (eye-index i i)) initial-element))
        eye)))
