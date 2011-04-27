;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla)

(defgeneric reconstruct (matrix-factorization)
  (:documentation "Calculate the original matrix from the matrix
  factorization/decomposition."))

(defclass lu ()
  ((lu :type matrix :initarg :lu :reader lu
       :documentation "matrix storing the transpose of the LU decomposition.")
   (ipiv :type vector :initarg :ipiv :reader ipiv
	 :documentation "pivot indices"))
  (:documentation "LU decomposition of a matrix with pivoting."))

(defclass qr ()
  ((qr :type matrix :initarg :qr :reader qr
       :documentation "matrix storing the QR decomposition."))
  (:documentation "QR decomposition of a matrix."))

(define-ondemand-slot ((qr qr) r)
  (bind (((:slots-r/o qr) qr)
         ((:accessors-r/o nrow ncol) qr))
    (assert (>= nrow ncol))
    (make-matrix :upper nil :initial-contents (matrix-from-first-rows qr ncol nil))))

;;; generic interface for square root-like decompositions

(defgeneric square-root (matrix-factorization left-or-right)
  (:documentation "Return the :LEFT or :RIGHT square root."))

;;; Cholesky decomposition

(defclass cholesky ()
  ((root :type (or lower-triangular-matrix upper-triangular-matrix)
         :initarg :root :reader root
         :documentation "Upper (lower) triangular matrix U (L) such that U^*U (LL^*)
                         is equal to the original matrix."))
  (:documentation "Cholesky decomposition a matrix."))

(defmethod initialize-instance :after ((instance cholesky) &key &allow-other-keys)
  (assert (typep (root instance) '(and triangular-matrix (satisfies square?)))))

(defmethod square-root ((cholesky cholesky) (left-or-right (eql :left)))
  (aetypecase (root cholesky)
    (lower-triangular-matrix it)
    (upper-triangular-matrix (transpose it))))

(defmethod square-root ((cholesky cholesky) (left-or-right (eql :right)))
  (aetypecase (root cholesky)
    (lower-triangular-matrix (transpose it))
    (upper-triangular-matrix it)))

(defmethod reconstruct ((cholesky cholesky))
  (aetypecase (root cholesky)
    (lower-triangular-matrix (mm it t))
    (upper-triangular-matrix (mm t it))))

;; (defgeneric permutations (object)
;;   (:documentation "Return the number of permutations in object (which is usually
;;   a matrix factorization, or a pivot index."))

;; (defun count-permutations% (ipiv)
;;   "Count the permutations in a pivoting vector."
;;   (iter
;;       (for index :from 1)               ; lapack counts from 1
;;       (for i :in-vector ipiv)
;;       (counting (/= index i))))

;; (defmethod permutations ((lu lu))
;;   (count-permutations% (ipiv lu)))



;; (defclass hermitian (matrix-factorization)
;;   ((factor :type (or lower-matrix upper-matrix)
;;            :initarg :factor :reader factor
;;            :documentation "upper/lower triangular matrix M such
;;              that MDM^* is equal to the original matrix")
;;    (ipiv :type vector :initarg :ipiv :reader ipiv
;;          :documentation "pivot indices"))
;;   (:documentation "Factorization for an indefinite hermitian matrix
;;   with pivoting."))

;; (defmethod component ((mf cholesky) component &key (copy? nil))
;;   (flet ((copy-maybe (matrix)
;;            (if copy?
;;                (copy-matrix matrix :copy? t)
;;                matrix)))
;;     (bind (((:slots-read-only factor) mf))
;;       (etypecase factor
;;         (lower-matrix
;;            (ecase component
;;              (:U (copy-maybe (conjugate-transpose factor)))
;;              (:L factor)))
;;         (upper-matrix
;;            (ecase component
;;              (:U factor)
;;              (:L (copy-maybe (conjugate-transpose factor)))))))))
