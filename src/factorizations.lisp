;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla)

(defclass lu ()
  ((lu :type matrix :initarg :lu :reader lu
       :documentation "matrix storing the transpose of the LU factorization.")
   (ipiv :type vector :initarg :ipiv :reader ipiv
	 :documentation "pivot indices"))
  (:documentation "LU factorization of a matrix with pivoting."))

(defun lu-u (lu)
  (make-upper-triangular-matrix (lu lu)))

(defun lu-l (lu)
  (aprog1 (make-lower-triangular-matrix (copy-array (lu lu)))
    (let+ (((&slots-r/o elements) it)
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
       :documentation "matrix storing the QR factorization.")
   (tau :accessor tau :initarg :tau :documentation "complex scalar for
   elementary reflectors (see documentation of xGEQRF)."))
  (:documentation "QR factorization of a matrix."))

(defun qr-r (qr &key copy?)
  (let+ (((&slots-r/o qr) qr)
         ((&accessors-r/o nrow ncol) qr))
    (assert (>= nrow ncol))
    (make-upper-triangular-matrix
     (clnu::maybe-copy-array (partition qr 0 ncol) copy?))))

;;; generic interface for square root-like factorizations

(defclass matrix-square-root ()
  ((left-square-root :reader left-square-root :initarg :left-square-root
                     :documentation "Matrix L such that LL^* is equal to the
                     original (decomposed) matrix.  This method should be
                     defined for other classes that can yield something
                     similar."))
  (:documentation "General class for representing all kinds of matrix square
  roots, regardless of how they were computed.  The convention is to store the
  left square root."))

(defun matrix-square-root (left-square-root)
  "Convenience function to create a matrix from a squre root."
  (make-instance 'matrix-square-root :left-square-root left-square-root))

(defmethod print-object ((matrix-square-root matrix-square-root) stream)
  (print-unreadable-object (matrix-square-root stream :type t)
    (format stream " LL^* with L=~A" (left-square-root matrix-square-root))))

(defgeneric right-square-root (object)
  (:documentation "Matrix L such that LL^* is equal to the
   original (decomposed) matrix.  Efficiency note: may be calculated on
   demand.")
  (:method (object)
    (transpose* (left-square-root object))))

(defmethod e2* ((a matrix-square-root) (b number))
  (make-instance (class-of a)
                 :left-square-root (e2* (left-square-root a) (sqrt b))))
(defmethod e2* ((a number) (b matrix-square-root))
  (e2* b a))
(defmethod e2/ ((a matrix-square-root) (b number))
  (make-instance (class-of a)
                 :left-square-root (e2/ (left-square-root a) (sqrt b))))


;;; Cholesky factorization

(defclass cholesky (matrix-square-root)
  ()
  (:documentation "Cholesky factorization a matrix."))

(defmethod initialize-instance :after ((instance cholesky) &key)
  (assert (typep (left-square-root instance) 
                 '(and lower-triangular-matrix (satisfies square?)))))

;;; permutations (pivoting)

(defgeneric permutations (object)
  (:documentation "Return the number of permutations in object (which is
  usually a matrix factorization, or a pivot index."))

(defun count-permutations% (ipiv)
  "Count the permutations in a pivoting vector."
  (iter
      (for index :from 1)               ; lapack counts from 1
      (for i :in-vector ipiv)
      (counting (/= index i))))

(defmethod permutations ((lu lu))
  (count-permutations% (ipiv lu)))

;;; hermitian factorization

(defclass hermitian-factorization ()
  ((factor :type matrix :initarg :factor :reader factor
           :documentation "see documentation of *SYTRF and *HETRF, storage is
           in the half specified by HERMITIAN-ORIENTATION and otherwise
           treated as opaque.")
   (ipiv :type vector :initarg :ipiv :reader ipiv :documentation "pivot
   indices"))
  (:documentation "Factorization for an indefinite hermitian matrix with
  pivoting."))

;;; spectral factorization

(defstruct spectral-factorization
  "Z W Z^T factorization of a Hermitian matrix, where the columns of Z contain
  the eigenvectors and W is a diagonal matrix of the eigenvalues.  Z is a
  unitary matrix."  z w)

;;; svd

(defstruct svd
  "Singular value decomposition.  Singular values are in S, in descending
order.  U and VT may be NIL in case they are not computed."
  (u nil) d (vt nil))

