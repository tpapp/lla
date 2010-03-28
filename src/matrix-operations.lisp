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
            :from (the fixnum (cm-index2 source-ld 0 col source-offset)))
       (for destination-index
            :from (cm-index2 destination-ld 0 col destination-offset))
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
dimensions and offsets.  If types are different, elements are
coerced."
  (if (eq source-type destination-type)
      (copy-columns% nrow ncol
                     source source-offset source-ld source-type
                     destination destination-offset destination-ld)
      (let ((destination-lisp-type (lla-type->lisp-type destination-type)))
        (copy-columns-loop% (coerce (row-major-aref source source-index)
                                    destination-lisp-type))))
  (values))

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

(defun stack-vertically (&rest arguments)
  "Stack arguments vertically, converting to a common type.  A vector
  is interpreted as a row matrix."
  (let* ((matrices (mapcar (lambda (arg)
                             (etypecase arg
                               (diagonal (diagonal->matrix arg))
                               (dense-matrix-like (set-restricted arg))
                               (numeric-vector (vector->row arg))))
                           arguments))
         (common-type (apply #'common-target-type (mapcar #'lla-type matrices)))
         (nrow-total (reduce #'+ matrices :key #'nrow))
         (ncol (ncol (first matrices)))
         (result (make-matrix common-type nrow-total ncol))
         (result-elements (elements result))
         (row 0))
    (iter
      (for matrix :in matrices)
      (for nrow := (nrow matrix))
      (for elements := (elements matrix))
      (for lla-type := (lla-type matrix))
      (assert (= ncol (ncol matrix)) ()
              "Matrix columns (or vector lengths) do not match.")
      (copy-columns nrow ncol 
                    elements 0 nrow lla-type
                    result-elements (cm-index2 nrow-total row 0) nrow-total common-type)
      (incf row nrow))
    result))

(defun stack-horizontally (&rest arguments)
  (declare (ignore arguments))
  (error "This function needs to be implemented."))
