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

(define-ondemand-slot ((lu lu) u)
  (make-matrix :upper nil :initial-contents (lu lu)))

(define-ondemand-slot ((lu lu) l)
  (aprog1 (make-matrix :lower nil :initial-contents (lu lu) :copy? t)
    (bind (((:slots-r/o elements) it)
           ((nrow ncol) (array-dimensions elements))
           (one (one* elements)))
      (dotimes (index (min nrow ncol))
        (setf (aref elements index index) one)))))

(defmethod print-object ((lu lu) stream)
  (print-unreadable-object (lu stream :type t)
    (with-slots (l u ipiv) lu
      (format stream "~2& L=~A~2& U=~A~2&  pivot indices=~A" l u ipiv))))

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

(defclass matrix-square-root ()
  ((left-square-root :reader left-square-root :initarg :left-square-root
                     :documentation "Matrix L such that LL^* is equal to the
                     original (decomposed) matrix.  This method should be defined for
                     other classes that can yield something similar."))
  (:documentation "General class for representing all kinds of matrix square roots,
  regardless of how they were computed.  The convention is to store the left square
  root."))

(defmethod print-object ((matrix-square-root matrix-square-root) stream)
  (print-unreadable-object (matrix-square-root stream :type t)
    (format stream " LL^* with L=~A" (left-square-root matrix-square-root))))

(defgeneric right-square-root (object)
  (:documentation "Matrix L such that LL^* is equal to the original (decomposed)
   matrix.  Efficiency note: may be calculated on demand.")
  (:method (object)
    (transpose* (left-square-root object))))

(defmethod reconstruct ((matrix-square-root matrix-square-root))
  (mm (left-square-root matrix-square-root) t))

(defmethod e2* ((a matrix-square-root) (b number))
  (make-instance (class-of a)
                 :left-square-root (e* (left-square-root a) (sqrt b))))
(defmethod e2* ((a number) (b matrix-square-root))
  (e2* b a))

;;; Cholesky decomposition

(defclass cholesky (matrix-square-root)
  ()
  (:documentation "Cholesky decomposition a matrix."))

(defmethod initialize-instance :after ((instance cholesky) &key &allow-other-keys)
  (assert (typep (left-square-root instance) 
                 '(and lower-triangular-matrix (satisfies square?)))))

;;; permutations (pivoting)

(defgeneric permutations (object)
  (:documentation "Return the number of permutations in object (which is usually
  a matrix factorization, or a pivot index."))

(defun count-permutations% (ipiv)
  "Count the permutations in a pivoting vector."
  (iter
      (for index :from 1)               ; lapack counts from 1
      (for i :in-vector ipiv)
      (counting (/= index i))))

(defmethod permutations ((lu lu))
  (count-permutations% (ipiv lu)))


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
