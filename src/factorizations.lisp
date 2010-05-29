;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla)

(define-abstract-class matrix-factorization ()
  ()
  (:documentation "Matrix factorization.  May not contain all
  components of the factorization."))

(defgeneric component (mf component &key copy?)
  (:documentation "Return a given component of a matrix factorization."))

(defgeneric reconstruct (mf)
  (:documentation "Calculate the original matrix from the matrix factorization."))

(defclass lu (matrix-factorization)
  ((lu-matrix :type dense-matrix :initarg :lu-matrix :reader lu-matrix
           :documentation "matrix storing the LU decomposition.")
   (ipiv :type vector :initarg :ipiv :reader ipiv
	 :documentation "pivot indices"))
  (:documentation "LU decomposition of a matrix with pivoting."))

(defclass qr (matrix-factorization)
  ((qr-matrix :type dense-matrix :initarg :qr-matrix :reader qr-matrix
           :documentation "matrix storing the QR decomposition."))
  (:documentation "QR decomposition of a matrix."))

(defclass cholesky (matrix-factorization)
  ((factor :type (or lower-matrix upper-matrix)
           :initarg :factor :reader factor
           :documentation "upper/lower triangular matrix U/L such
             that U^*U or LL^* is equal to the original matrix"))
  (:documentation "Cholesky decomposition a matrix."))

(defmethod initialize-instance :after ((instance cholesky) &key &allow-other-keys)
  (assert (typep (factor instance) '(and square-matrix
                                     (or lower-matrix
                                      upper-matrix)))))

(defclass hermitian (matrix-factorization)
  ((factor :type (or lower-matrix upper-matrix)
           :initarg :factor :reader factor
           :documentation "upper/lower triangular matrix M such
             that MDM^* is equal to the original matrix")
   (ipiv :type vector :initarg :ipiv :reader ipiv
         :documentation "pivot indices"))
  (:documentation "Factorization for an indefinite hermitian matrix
  with pivoting."))

(defmethod component ((mf qr) (component (eql :R)) &key copy?)
  (declare (ignore copy?))
  (bind (((:slots-read-only nrow ncol elements) (qr-matrix mf)))
    (matrix-from-first-rows elements ncol ncol nrow :upper)))

(defmethod component ((mf cholesky) component &key (copy? nil))
  (flet ((copy-maybe (matrix)
           (if copy?
               (copy-matrix matrix :copy? t)
               matrix)))
    (bind (((:slots-read-only factor) mf))
      (etypecase factor
        (lower-matrix
           (ecase component
             (:U (copy-maybe (conjugate-transpose factor)))
             (:L factor)))
        (upper-matrix
           (ecase component
             (:U factor)
             (:L (copy-maybe (conjugate-transpose factor)))))))))
