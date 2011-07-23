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

(defstruct (matrix-square-root (:constructor make-matrix-square-root (left)))
  "General class for representing XX^T decompositions of matrices, regardless
  of how they were computed.  The convention is to store X, the left square
  root."
  left)

(defgeneric left-square-root (a)
  (:documentation "Return X such that XX^T=A.")
  (:method ((a matrix-square-root))
    (matrix-square-root-left a)))

(defgeneric right-square-root (a)
  (:documentation "Return Y such that Y^T Y=A.  Efficiency note:
  decompositions should store the left square root X, and compute Y=X^T on
  demand, so getting X directly might be more efficient.")
  (:method ((a matrix-square-root))
    (transpose (matrix-square-root-left a))))

(declaim (inline xx))
(defun xx (left-square-root)
  "Convenience function to create a matrix from a left square root."
  (make-matrix-square-root left-square-root))

(defmacro define-matrix-square-root-scalar-multiplication
    (type &key (make (symbolicate '#:make- type)))
  `(progn
     (defmethod e2* ((a ,type) (b number))
       (,make (e2* (left-square-root a) (sqrt b))))
     (defmethod e2* ((a number) (b ,type))
       (,make (e2* (sqrt a) (left-square-root b))))
     (defmethod e2/ ((a ,type) (b number))
       (,make (e2/ (left-square-root a) (sqrt b))))))

(define-matrix-square-root-scalar-multiplication matrix-square-root)

;;; Cholesky factorization

(defstruct (cholesky (:include matrix-square-root)
                     (:constructor make-cholesky% (left)))
  "Cholesky factorization a matrix.")

(defun make-cholesky (left)
  (assert (typep left '(and lower-triangular-matrix (satisfies square?))))
  (make-cholesky% left))

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
  unitary matrix." z w)

;;; svd

(defstruct svd
  "Singular value decomposition.  Singular values are in S, in descending
order.  U and VT may be NIL in case they are not computed."
  (u nil) d (vt nil))
