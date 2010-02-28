(in-package :lla)

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
                    (make-matrix* ,lla-type ncol nrow transposed :kind transposed-matrix-kind))))))

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
         (copy-matrix matrix :copy-p t)))))

;;;; extraction methods for factorization components and related
;;;; utility functions

(defun copy-columns% (nrow ncol source source-type ld-source
                      destination destination-type ld-destination
                      &optional (destination-offset 0))
  "Copy columns from source to destination."
  (dotimes (col ncol)
    (copy-elements-into source source-type (* col ld-source)
                        destination destination-type (+ destination-offset (* col ld-destination))
                        nrow))
  (values))

(defun matrix-from-first-rows (lla-type vector nrow ncol leading-dimension &optional (kind :dense))
  "Extract & return (as KIND, default is :DENSE) the first NROW rows
of an LEADING-DIMENSIONxNCOL matrix, given as a Lisp vector with
column-major indexing.  NOTE: needed to interface to LAPACK routines
like xGELS."
  (let ((result (make-matrix lla-type nrow ncol :kind kind)))
    (copy-columns% nrow ncol vector lla-type leading-dimension
                   (elements result) lla-type nrow)
    result))

(defmethod component ((mf qr) (component (eql :R)) &key copy-p)
  (declare (ignore copy-p))
  (bind (((:slots-read-only qr-matrix) mf)
         ((:slots-read-only nrow ncol elements) qr-matrix))
    (matrix-from-first-rows (lla-type qr-matrix) elements ncol ncol nrow :upper-triangular)))

(defmethod component ((mf cholesky) component &key (copy-p nil))
  (flet ((copy-maybe (matrix)
           (if copy-p
               (copy-matrix matrix :copy-p t)
               matrix)))
    (bind (((:slots-read-only factor) mf))
      (etypecase factor
        (lower-triangular-matrix
           (ecase component
             (:U (copy-maybe (transpose factor)))
             (:L factor)))
        (upper-triangular-matrix
           (ecase component
             (:U factor)
             (:L (copy-maybe (transpose factor)))))))))

     
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
      (assert (= ncol (ncol matrix)) () "Matrix columns (or vector lengths) do not match.")
      (iter
        (for col :from 0 :below ncol)
        (copy-elements-into elements lla-type (cm-index2 nrow 0 col)
                            result-elements common-type (cm-index2 nrow-total row col)
                            nrow))
      (incf row nrow))
    result))

(defun stack-horizontally (&rest arguments)
  (declare (ignore arguments))
  (error "This function needs to be implemented."))
