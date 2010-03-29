(in-package :lla)

;;;  Matrix classes --- abstract interface
;;;
;;;  Being a subclass of DENSE-MATRIX-LIKE essentially specifies a
;;;  storage model: elements are stored in column major format, and
;;;  the resulting vector can be passed to LAPACK.  See SET-RESTRICTED
;;;  and RESTRICTED-ELEMENTS for nuances.

(define-abstract-class dense-matrix-slots ()
  ((nrow :type dimension :initarg :nrow :reader nrow
	 :documentation "The number of rows in the matrix.")
   (ncol :type dimension :initarg :ncol :reader ncol
	 :documentation "The number of columns in the matrix."))
  (:documentation "Defines NROW & NCOL slots and associated XDIM,
  XDIMS and XSIZE queries.  Meant to be used as a mixin class, implies
  no element mapping."))

(defmethod xrank ((matrix dense-matrix-slots))
  2)

(defmethod xdims ((matrix dense-matrix-slots))
  (list (nrow matrix) (ncol matrix)))

(defmethod xdim ((matrix dense-matrix-slots) axis-number)
  (ecase axis-number
    (0 (nrow matrix))
    (1 (ncol matrix))))

(defmethod xsize ((matrix dense-matrix-slots))
  (* (nrow matrix) (ncol matrix)))

(define-abstract-class dense-matrix-like (dense-matrix-slots numeric-vector-like)
  ()
  (:documentation "Superclass of all classes that resemble a matrix.
  Nevertheless, it is not implied that the elements should be
  interpreted as a matrix: they could be a matrix factorization packed
  into a matrix, etc.

  Subclasses of this class behave like an NROW x NCOL matrix, storing
  the elements in a subset of ELEMENTS, mapped using _column-major_
  indexing.  Specifically, the position of a particular element at
  index (row,col) row + NCOL*col, with 0 <= row < NROW and 0 <= col <
  NCOL.  The function CM-INDEX2 can be used for calculations.

  If the subclass is also a member of RESTRICTED-ELEMENTS, then not
  all elements are necessarily stored: some may be constants (eg for
  trianguar matrices), some may be inferred from other elements (eg
  for hermitian matrices.  When SET-RESTRICTED is called, it has to
  ensure that all values in ELEMENTS reflect the constant and/or
  inferred cells in the matrix."))

(declaim (inline cm-index2)
         (ftype (function (fixnum fixnum fixnum) fixnum)))
(defun cm-index2 (nrow row col)
  "Calculate column-major index, without error checking.  Inlined."
  (the fixnum (+ (the fixnum (* nrow col)) row)))

;;; General XREF and (SETF XREF) methods for all DENSE-MATRIX-LIKE
;;; objects.

(defmethod xref ((matrix dense-matrix-like) &rest subscripts)
  (bind (((row col) subscripts)
         ((:accessors-r/o elements nrow ncol) matrix))
    (check-index row nrow)
    (check-index col ncol)
    (aref elements (cm-index2 nrow row col))))

(defmethod (setf xref) (value (matrix dense-matrix-like) &rest subscripts)
  (bind (((row col) subscripts)
         ((:accessors-r/o elements nrow ncol) matrix))
    (check-index row nrow)
    (check-index col ncol)
    (setf (aref elements (cm-index2 nrow row col))
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
