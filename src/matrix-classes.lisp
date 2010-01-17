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
  (:documentation "Superclass of all classes that resemble a matrix.
  Nevertheless, it is not implied that the elements should be
  interpreted as a matrix: they could be a matrix factorization packed
  into a matrix, etc.  Subclasses of this class behave like an NROW x
  NCOL matrix, storing the elements in a subset of ELEMENTS, mapped
  using column-major indexing.  If the subclass is also a member of
  RESTRICTED-ELEMENTS, then not all elements are necessarily stored:
  some may be constants (eg for trianguar matrices), some may be
  inferred from other elements (eg for hermitian matrices.  When
  SET-RESTRICTED is called, it has to ensure that all values in
  ELEMENTS reflect the constant or inferred cells in the matrix."))

(define-abstract-class dense-matrix-like-square (dense-matrix-like)
  ()
  (:documentation "Enforces square dimensions."))

(defmethod initialize-instance :after ((object dense-matrix-like-square)
                                       &key &allow-other-keys)
  (assert (= (nrow object) (ncol object))))

(defmacro define-matrix-class ()
  "Define the function MATRIX-CLASS."
  (let ((kind-prefix-pairs '((:dense dense-matrix)
                             (:lower-triangular lower-triangular-matrix)
                             (:upper-triangular upper-triangular-matrix)
                             (:hermitian hermitian-matrix))))
    (labels ((generate-lla-type-case (prefix)
               `(ecase lla-type
                  ,@(mapcar (lambda (lla-type)
                              `(,lla-type ',(make-symbol* prefix '- lla-type)))
                       +lla-type-list+)))
             (generate-kind-case ()
               `(ecase kind
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

(defmacro define-dense-matrix-subclass
    (kind (&rest additional-superclasses) documentation
     &optional (additional-slots '()))
  "Macro for defining simple subclasses of DENSE-MATRIX-LIKE."
  (check-type kind symbol)
  (check-type documentation string)
  (let ((class (make-symbol* kind '-matrix)))
    `(progn
       (define-abstract-class ,class (,@additional-superclasses dense-matrix-like)
         ,additional-slots
         (:documentation ,documentation))
       (define-lla-class ,class)
       (defmethod matrix-kind ((matrix ,class)) ,kind))))

(define-dense-matrix-subclass :dense ()
  "Dense matrix, with elements stored in column-major order.")

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
  (:method ((matrix numeric-vector))
    ;; do nothing, return matrix
    matrix)
  (:documentation "Set restricted/unaccessed elements to the
appropriate value in the data vector of matrix.  Always return the
matrix.  Useful when calling functions which expect a proper dense
matrix."))

(define-dense-matrix-subclass :upper-triangular (restricted-elements)
    "A dense, upper triangular matrix.  The elements below the
diagonal are not necessarily initialized and not accessed.")

(define-dense-matrix-subclass :lower-triangular (restricted-elements)
    "A dense, lower triangular matrix.  The elements above the
diagonal are not necessarily initialized and not accessed.")

(define-dense-matrix-subclass :hermitian
    (restricted-elements dense-matrix-like-square)
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
  (:documentation "Matrix factorization.  May not contain all
  components of the factorization."))

(defgeneric factorization-component (mf component &key copy-p)
  (:documentation "Return a given component of a matrix factorization."))

(defgeneric reconstruct (mf)
  (:documentation "Calculate the original matrix from the matrix factorization."))

(defclass lu (matrix-factorization)
  ((lu-matrix :type dense-matrix :initarg :lu-matrix :reader lu-matrix
           :documentation "matrix storing the LU decomposition.")
   (ipiv :type numeric-vector-integer :initarg :ipiv :reader ipiv
	 :documentation "pivot indices"))
  (:documentation "LU decomposition of a matrix with pivoting."))

(defclass qr (matrix-factorization)
  ((qr-matrix :type dense-matrix :initarg :qr-matrix :reader qr-matrix
           :documentation "matrix storing the QR decomposition."))
  (:documentation "QR decomposition of a matrix."))

(defclass cholesky (matrix-factorization)
  ((r-matrix :type upper-triangular-matrix :initarg :r-matrix :reader r-matrix
             :documentation "upper triangular matrix R such that R^*R
             is equal to the original matrix"))
  (:documentation "Cholesky decomposition a matrix."))

(defmethod initialize-instance :after ((instance cholesky) &key &allow-other-keys)
  (assert (typep (r-matrix instance) 'upper-triangular-matrix)))
