(in-package :lla)

;;;; Matrix classes --- abstract interface
;;;
;;; We define "abstract" superclasses which are then specialized
;;; according to LLA-TYPE.  Each class inherits from its abstract
;;; superclass (identified with the KIND of matrix), and also from
;;; NUMERIC-VECTOR-LLA-TYPE.  DEFINE-LLA-CLASS is a helper macro for
;;; this.
;;;
;;; Being a subclass of DENSE-MATRIX-LIKE essentially specifies a
;;; storage model: elements are stored in column major format, and the
;;; resulting vector can be passed to LAPACK.  See SET-RESTRICTED and
;;; RESTRICTED-ELEMENTS for nuances.

(define-abstract-class dense-matrix-like (numeric-vector-like)
  ((nrow :type dimension :initarg :nrow :reader nrow
	 :documentation "The number of rows in the matrix.")
   (ncol :type dimension :initarg :ncol :reader ncol
	 :documentation "The number of columns in the matrix."))
  (:documentation "Superclass of all classes that resemble a matrix.
  Nevertheless, it is not implied that the elements should be
  interpreted as a matrix: they could be a matrix factorization packed
  into a matrix, etc.

  Subclasses of this class behave like an NROW x NCOL matrix, storing
  the elements in a subset of ELEMENTS, mapped using column-major
  indexing.  If the subclass is also a member of RESTRICTED-ELEMENTS,
  then not all elements are necessarily stored: some may be
  constants (eg for trianguar matrices), some may be inferred from
  other elements (eg for hermitian matrices.  When SET-RESTRICTED is
  called, it has to ensure that all values in ELEMENTS reflect the
  constant or inferred cells in the matrix.

  The length of ELEMENTS does not have to equal (* NROW NCOL), but
  (* (LEADING-DIMENSION MATRIX) NCOL), this is to allow for adjustable
  matrices and at the moment it is not used for anything else (such as
  views, etc)."))

(defgeneric leading-dimension (matrix)
  (:documentation "The leading dimension of the matrix.")
  (:method ((matrix dense-matrix-like))
    (nrow matrix)))

(defun square-matrix-p (matrix)
  "Test if a matrix is square."
  (= (nrow matrix) (ncol matrix)))

(deftype square-matrix ()
  '(and dense-matrix-like
    (satisfies square-matrix-p)))

(defgeneric matrix-class (kind)
  (:documentation "Return the name of the matrix class corresponding
  to KIND (which can be :dense, :upper-triangular, etc)."))

(defgeneric matrix-kind (matrix)
  (:documentation "Return the matrix kind, eg :DENSE, :UPPER-TRIANGULAR, etc."))

(defmacro define-dense-matrix-subclass
    (kind (&rest additional-superclasses) documentation
     &optional (additional-slots '()))
  "Macro for defining simple subclasses of DENSE-MATRIX-LIKE."
  (check-type kind symbol)
  (check-type documentation string)
  (let ((class (make-symbol* kind '-matrix))
        (kind (make-keyword* kind)))
    `(progn
       (defclass ,class (,@additional-superclasses dense-matrix-like)
         ,additional-slots
         (:documentation ,documentation))
       (defmethod matrix-kind ((matrix ,class)) ,kind)
       (defmethod matrix-class ((kind (eql ,kind))) ',class))))

;;; Framework for dense matrices with restricted elements.

(define-abstract-class restricted-elements ()
  ()
  (:documentation "Superclass for objects with restricted elements.  A
  restricted element is not used in ELEMENTS (eg it is either a
  constant like 0, or calculated for other elements).  Nevertheless,
  it needs to be set when the matrix is passed to certain LAPACK
  routines.  Calling SET-RESTRICTED will ensure that it is set to the
  correct value."))

(defgeneric set-restricted (matrix)
  (:method ((matrix numeric-vector-like))
    ;; do nothing, return matrix
    matrix)
  (:documentation "Set restricted/unaccessed elements to the
appropriate value in the data vector of matrix.  Always return the
matrix.  Useful when calling functions which expect a proper dense
matrix."))

(define-dense-matrix-subclass dense ()
  "Dense matrix, with elements stored in column-major order.")

(define-dense-matrix-subclass upper-triangular (restricted-elements)
    "A dense, upper triangular matrix.  The elements below the
diagonal are not necessarily initialized and not accessed.")

(define-dense-matrix-subclass lower-triangular (restricted-elements)
    "A dense, lower triangular matrix.  The elements above the
diagonal are not necessarily initialized and not accessed.")

(define-dense-matrix-subclass hermitian (restricted-elements)
  ;; LLA uses the class HERMITIAN-MATRIX to implement both real
  ;; symmetric and complex Hermitian matrices --- as technically, real
  ;; symmetric matrices are also Hermitian.  Complex symmetric
  ;; matrices are NOT implemented as a special matrix type, as they
  ;; don't have any special properties (eg real eigenvalues, etc).
  "A dense Hermitian matrix, with elements stored in the upper
  triangle.")

(defmethod initialize-instance :after ((object hermitian-matrix)
                                       &key &allow-other-keys)
  (check-type object square-matrix))


;;; Matrix factorizations

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
