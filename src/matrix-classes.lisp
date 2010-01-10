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

(define-abstract-class dense-matrix-like ()
  ((nrow :type dimension :initarg :nrow :reader nrow
	 :documentation "The number of rows in the matrix.")
   (ncol :type dimension :initarg :ncol :reader ncol
	 :documentation "The number of columns in the matrix."))
  (:documentation "Subclasses of this class behave like an NROW x NCOL
  matrix, storing the elements in a subset of ELEMENTS, mapped using
  column-major indexing.  If the subclass is also a member of
  RESTRICTED-ELEMENTS, then not all elements are necessarily stored:
  some may be constants (eg for trianguar matrices), some may be
  inferred from other elements (eg for hermitian matrices.  When
  SET-RESTRICTED is called, it has to ensure that all values in
  ELEMENTS reflect the constant or inferred cells in the matrix."))

(defmacro define-matrix-class ()
  "Define the function MATRIX-CLASS."
  (let ((kind-prefix-pairs '((:dense dense-matrix)
                             (:lower-triangular lower-triangular-matrix)
                             (:upper-triangular upper-triangular-matrix)
                             (:hermitian hermitian-matrix)
                             (:lu lu-factorization)
                             (:qr qr-factorization)
                             (:cholesky cholesky-factorization))))
    (labels ((generate-lla-type-case (prefix)
               `(case lla-type
                  ,@(mapcar (lambda (lla-type)
                              `(,lla-type ',(make-symbol* prefix '- lla-type)))
                       +lla-type-list+)))
             (generate-kind-case ()
               `(case kind
                  ,@(mapcar (lambda (kind-and-prefix)
                              `(,(first kind-and-prefix)
                                 ,(generate-lla-type-case (second kind-and-prefix))))
                       kind-prefix-pairs))))
      `(defun matrix-class (kind lla-type)
         "Return symbol for matrix class."
         ,(generate-kind-case)))))

(declaim (inline matrix-class))
(define-matrix-class)


(defgeneric matrix-kind (matrix)
  (:documentation "Return the matrix kind,
  eg :DENSE, :UPPER-TRIANGULAR, etc."))

(defmacro define-dense-matrix-kind ((kind &optional (class-name (make-symbol* kind '-matrix)))
                                    (&rest additional-superclasses) documentation
                                    &optional (additional-slots '()))
  "Macro for defining simple subclasses of DENSE-MATRIX-LIKE."
  (check-type kind symbol)
  (check-type class-name symbol)
  (check-type documentation string)
  `(progn
     (define-abstract-class ,class-name (,@additional-superclasses dense-matrix-like)
       ,additional-slots
       (:documentation ,documentation))
     (define-lla-class ,class-name)
     (defmethod matrix-kind ((matrix ,class-name)) ,kind)))

(define-dense-matrix-kind (:dense) ()
  "Dense matrix, with elements stored in column-major order.")

;;; Framework for dense matrices with restricted elements.

(define-abstract-class restricted-elements ()
  ((restricted-set-p :type boolean :accessor restricted-set-p
                     :initarg :restricted-set-p
                     :initform nil
                     :documentation "non-nil iff not accessible
   elements in the data vectors have been enforced to contain 0."))
  (:documentation "Superclass for objects with restricted elements.  A
  restricted element is not used in ELEMENTS (eg it is either a
  constant like 0, or calculated for other elements).  Nevertheless,
  it needs to be set when the matrix is passed to certain LAPACK
  routines.  Calling SET-RESTRICTED will ensure that it is set to the
  correct value."))

(defgeneric set-restricted* (matrix)
  (:documentation "Always set the restricted elements of matrix,
  regardless of restricted-set-p.  Should not be called directly,
  use set-restricted instead."))

(defgeneric set-restricted (matrix)
  (:documentation "Set restricted/unaccessed elements to the
  appropriate value in the data vector of matrix.  Useful when calling
  functions which expect a proper dense matrix.  If restricted-set-p,
  or if the object is not a subtype of restricted-elements, do
  nothing.

  NOTE: the default behavior is to call set-restricted* if (not
  restricted-set-p), so it is advised to define that method instead
  for classes.")
  (:method ((matrix dense-matrix))) ;; do nothing
  (:method ((matrix restricted-elements))
    (with-slots (restricted-set-p) matrix
      (unless restricted-set-p
        (set-restricted* matrix)
        (setf restricted-set-p t)))))

(defmethod restricted-set-p ((matrix dense-matrix))
  ;; Dense matrices can behave as if their restricted elements were
  ;; always set.  This makes handling some special cases easier.
  t)

(define-dense-matrix-kind (:upper-triangular) (restricted-elements)
    "A dense, upper triangular matrix.  The elements below the
diagonal are not necessarily initialized and not accessed.")

(define-dense-matrix-kind (:lower-triangular) (restricted-elements)
    "A dense, lower triangular matrix.  The elements above the
diagonal are not necessarily initialized and not accessed.")

(define-dense-matrix-kind (:hermitian) (restricted-elements)
  ;; LLA uses the class HERMITIAN-MATRIX to implement both real
  ;; symmetric and complex Hermitian matrices --- as technically, real
  ;; symmetric matrices are also Hermitian.  Complex symmetric
  ;; matrices are NOT implemented as a special matrix type, as they
  ;; don't have any special properties (eg real eigenvalues, etc).
  "A dense Hermitian matrix, with elements stored in the upper
  triangle.")


;;; Matrix factorizations

(define-abstract-class matrix-factorization ()
  ()
  (:documentation "Technically these are not matrices (as they may
contain 2+ matrices and other data), but they are usually defined as
subclasses of dense-matrix-like (or a subclass, etc) so that you can
manipulate them as if they were.  NROW and NCOL follow the conventions
of LAPACK, ie usually refer to the dimensions of the original
matrix."))

(defgeneric factorization-component (mf component)
  (:documentation "Return a given component of a matrix factorization."))

(defgeneric reconstruct (mf)
  (:documentation "Calculate the original matrix from the matrix factorization."))

(define-dense-matrix-kind (:lu lu-factorization) (matrix-factorization)
  "LU decomposition of the matrix A with pivoting."
  ((ipiv :type numeric-vector-integer :initarg :ipiv :reader ipiv
	 :documentation "pivot indices")))

(define-dense-matrix-kind (:qr qr-factorization) (matrix-factorization)
  "QR decomposition of a matrix.")

(define-dense-matrix-kind (:cholesky
cholesky-factorization) (matrix-factorization lower-triangular-matrix)
  "Cholesky decomposition R of a matrix=RR^*, where * is the conjugate
transpose.  Should behave as a matrix, but only the lower-triangular
part of the decomposition is stored.")

