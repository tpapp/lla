(in-package :lla)

;;;  Matrix classes --- abstract interface
;;;
;;;  We define "abstract" superclasses which are then specialized
;;;  according to LLA-TYPE.  Each class inherits from its abstract
;;;  superclass (identified with the KIND of matrix), and also from
;;;  NUMERIC-VECTOR-LLA-TYPE.  DEFINE-LLA-CLASS is a helper macro for
;;;  this.
;;;
;;;  Being a subclass of DENSE-MATRIX-LIKE essentially specifies a
;;;  storage model: elements are stored in column major format, and
;;;  the resulting vector can be passed to LAPACK.  See SET-RESTRICTED
;;;  and RESTRICTED-ELEMENTS for nuances.

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
  the elements in a subset of ELEMENTS, mapped using _column-major_
  indexing.  Specifically, the position of a particular element at
  index (row,col) is OFFSET + row + LEADING-DIMENSION*col, with 0 <=
  row < NROW and 0 <= col < NCOL.  The function CM-INDEX2 can be used
  for calculations.  OFFSET and LEADING-DIMENSION are not necessarily
  slots in a particular (sub)class, they should be queried using
  accessors.

  If the subclass is also a member of RESTRICTED-ELEMENTS, then not
  all elements are necessarily stored: some may be constants (eg for
  trianguar matrices), some may be inferred from other elements (eg
  for hermitian matrices.  When SET-RESTRICTED is called, it has to
  ensure that all values in ELEMENTS reflect the constant and/or
  inferred cells in the matrix."))

(define-abstract-class compact-matrix ()
  ()
  (:documentation "Compact storage model for matrices, with
  LEADING-DIMENSION=NROW and OFFSET=0"))

(defgeneric leading-dimension (matrix)
  (:documentation "The leading dimension of the matrix.")
  (:method ((matrix compact-matrix))
    (nrow matrix)))

(defgeneric offset (matrix)
  (:documentation "Offset of the first element.")
  (:method ((matrix compact-matrix))
    0))

(declaim (ftype (function (t) dimension) leading-dimension offset))

(declaim (inline cm-index2)
         (ftype (function (fixnum fixnum fixnum &optional fixnum) fixnum)))
(defun cm-index2 (leading-dimension row col &optional (offset 0))
  "Calculate column-major index, without error checking.  Inlined."
  (the fixnum (+ (the fixnum (* leading-dimension col)) row offset)))

;;; General XREF and (SETF XREF) methods for all DENSE-MATRIX-LIKE
;;; objects.

(defmethod xref ((matrix dense-matrix-like) &rest subscripts)
  (bind (((row col) subscripts)
         ((:accessors-r/o elements nrow ncol leading-dimension offset) matrix))
    (check-index row nrow)
    (check-index col ncol)
    (aref elements (cm-index2 leading-dimension row col offset))))

(defmethod (setf xref) (value (matrix dense-matrix-like) &rest subscripts)
  (bind (((row col) subscripts)
         ((:accessors-r/o elements nrow ncol leading-dimension offset) matrix))
    (check-index row nrow)
    (check-index col ncol)
    (setf (aref elements (cm-index2 leading-dimension row col offset))
          value)))

;;;  Square matrix type

(defun square-matrix-p (matrix)
  "Test if a matrix is square."
  (= (nrow matrix) (ncol matrix)))

(deftype square-matrix ()
  '(and dense-matrix-like
    (satisfies square-matrix-p)))

;;;  Matrix kind and class
;;;
;;;  Matrices have a _kind_, eg dense, lower-trianguar, etc, with a
;;;  bijection between kinds and class names.  The following two
;;;  functions provide the framework for the mapping, and the macro
;;;  defines them automatically.

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
