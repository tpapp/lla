;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla)

(define-abstract-class matrix-factorization ()
  ()
  (:documentation "Matrix factorization.  May not contain all
  components of the factorization."))

(defgeneric component (mf component &key copy-p)
  (:documentation "Return a given component of a matrix factorization."))

(defgeneric reconstruct (mf)
  (:documentation "Calculate the original matrix from the matrix factorization."))

(defclass lu (matrix-factorization)
  ((lu-matrix :type dense-matrix :initarg :lu-matrix :reader lu-matrix
           :documentation "matrix storing the LU decomposition.")
   (ipiv :type numeric-vector :initarg :ipiv :reader ipiv
	 :documentation "pivot indices"))
  (:documentation "LU decomposition of a matrix with pivoting."))

(defclass qr (matrix-factorization)
  ((qr-matrix :type dense-matrix :initarg :qr-matrix :reader qr-matrix
           :documentation "matrix storing the QR decomposition."))
  (:documentation "QR decomposition of a matrix."))

(defclass cholesky (matrix-factorization)
  ((factor :type (or lower-triangular-matrix upper-triangular-matrix)
           :initarg :factor :reader factor
             :documentation "upper/lower triangular matrix U/L such
             that U^*U or LL^* is equal to the original matrix"))
  (:documentation "Cholesky decomposition a matrix."))

(defmethod initialize-instance :after ((instance cholesky) &key &allow-other-keys)
  (assert (typep (factor instance) '(and square-matrix
                                     (or lower-triangular-matrix
                                         upper-triangular-matrix)))))

(defclass hermitian-factorization (matrix-factorization)
  ((factor :type (or lower-triangular-matrix upper-triangular-matrix)
           :initarg :factor :reader factor
           :documentation "upper/lower triangular matrix M such
             that MDM^* is equal to the original matrix")
   (ipiv :type numeric-vector :initarg :ipiv :reader ipiv
         :documentation "pivot indices"))
  (:documentation "Factorization for an indefinite hermitian matrix
  with pivoting."))

(defmethod component ((mf qr) (component (eql :R)) &key copy-p)
  (declare (ignore copy-p))
  (bind (((:slots-read-only qr-matrix) mf)
         ((:slots-read-only nrow ncol elements) qr-matrix))
    (matrix-from-first-rows (lla-type qr-matrix) elements ncol ncol nrow :upper-triangular)))

(defmethod component ((mf cholesky) component &key (copy-p nil))
  (flet ((copy-maybe (matrix)
           (if copy-p
               (copy-matrix% matrix :copy-p t)
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
