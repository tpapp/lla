(in-package :lla)

;;;;  ??? currently, the interface is evolving in the direction of a
;;;;  dual semantics for factorizations and other special types that
;;;;  carry extra information.  This means that, for example, an LU
;;;;  decomposition is at the same time a matrix-like object (you can
;;;;  access elements as usual, etc), and something that carries extra
;;;;  information and behaves differently in extra information: (solve
;;;;  a b) would either return A^-1 B, or (LUP)^-1 B, even though LU
;;;;  might "look" the same when accessed with eg. xref. -- Tamas
;;;;
;;;;  !!! provide a take interface for conversions.

(define-abstract-class matrix ()
  ((nrow :type dimension :initarg :nrow :reader nrow
	 :documentation "The number of rows in the matrix.")
   (ncol :type dimension :initarg :ncol :reader ncol
	 :documentation "The number of columns in the matrix.")
   (data :type numeric-vector :initarg :data
	 :reader data
	 :documentation "Slot for the elements of the matrix.  Type is
	 always numeric-vector, but the length and interpretation of
	 contents is implementation-dependent.")))

;;; copy data on demand
(defmethod copy-data ((matrix matrix))
  (copy-data (data matrix)))

;;; LLA element type, from the type of data
(defmethod lla-type ((matrix matrix))
  (lla-type (data matrix)))

;;; Some of the xarray interface, xref defined later for specific
;;; subclasses

(defmethod xelttype ((matrix matrix))
  (xelttype (data matrix)))

(defmethod xrank ((matrix matrix))
  2)

(defmethod xdims ((matrix matrix))
  (list (nrow matrix) (ncol matrix)))

(defmethod xdim ((matrix matrix) axis-number)
  (ecase axis-number
    (0 (nrow matrix))
    (1 (ncol matrix))))

(defmethod xsize ((matrix matrix))
  (* (nrow matrix) (ncol matrix)))

;;; printing

(defvar *print-matrix-aligned* t
  "If non-nil, characters will be aligned.")

(defvar *print-matrix-paddig* 1
  "Number of spaces between columns.")

(defun upper-triangular-mask (row col)
  (if (<= row col)
      nil
      "."))

(defun lower-triangular-mask (row col)
  (if (>= row col)
      nil
      "."))

(defun print-matrix (matrix stream &key 
		     (formatter #'standard-numeric-formatter)
                     (accessor #'xref)
                     (mask (constantly nil)))
  "Format and print the elements of matrix to stream, using
*print-matrix-padding* spaces between columns.  If
*print-matrix-aligned*, columns will be right-aligned.  Prints at most
*print-length* rows and columns, indicating more with a ...  Mask is a
function called with row and col indexes, and should return nil for
elements printed.  If mask returns a non-nil value, that will be
printed instead (should be a string)."
  ;; ?? maybe column & row labels, not a high priority at the moment
  (bind (((:values nrow row-trunc-p) (print-length-truncate (nrow matrix)))
	 ((:values ncol col-trunc-p) (print-length-truncate (ncol matrix)))
	 (formatted-elements (make-array (list nrow ncol)))
	 (column-widths (make-array ncol :element-type 'fixnum :initial-element 0))
	 (padding (make-array *print-matrix-paddig*
                              :element-type 'character :initial-element #\space))
	 (aligned-p *print-matrix-aligned*))
    ;; first pass - format elements, measure width
    (dotimes (col ncol)
      (dotimes (row nrow)
	(let* ((mask (funcall mask row col))
               (formatted-element (if mask
                                      mask
                                      (funcall formatter 
                                               (funcall accessor matrix row col))))
	       (width (length formatted-element)))
	  (setf (aref column-widths col) (max (aref column-widths col) width)
		(aref formatted-elements row col) formatted-element))))
    ;; second pass - print
    (dotimes (row nrow)
      (when (plusp row)
	(fresh-line stream))
      (dotimes (col ncol)
	(when (plusp col)
	  (princ padding stream))
	(let ((elt (aref formatted-elements row col)))
	  (if aligned-p
	      (format stream "~V@A" (aref column-widths col) elt)
	      (princ elt stream))))
      (when col-trunc-p
	(princ padding stream)
	(princ "..." stream)))
    (when row-trunc-p
      (format stream "~&..."))))

(defmacro define-matrix-like-printer ((instance class stream) slots &body body)
  "Simple wrapper for the standard printing style."
  (check-type class symbol)
  (check-type instance symbol)
  (check-type stream symbol)
  `(defmethod print-object ((,instance ,class) ,stream)
     (print-unreadable-object (,instance ,stream :type t)
       (format ,stream "~a x ~a ~a elements~&" (nrow ,instance) (ncol ,instance)
               (lla-type ,instance))
       (with-slots ,slots ,instance
         ,@body))))

(define-matrix-like-printer (matrix matrix stream) ()
  (print-matrix matrix stream))

;;;;
;;;;  Dense matrix class
;;;;
;;;;  Dense matrices store elements in a column-major order.  Some
;;;;  other matrix classes also employ this storage scheme, but may
;;;;  not use all elements (eg triangular matrices).  The common
;;;;  superclass of all of these is dense-matrix-like.
;;;;
;;;;  IMPORTANT: whenever an argument is dense-matrix-like but need to
;;;;  be interpreted in a FORTRAN call as a dense-matrix, element
;;;;  restrictions need to be enforced.  The best way to do this is
;;;;  using (set-restricted object): this will ensure that next time
;;;;  this is done, it will not result in efficiency loss (unless the
;;;;  object is modified, of course).

(define-abstract-class dense-matrix-like (matrix)
  ()
  (:documentation "Elements stored like a dense matrix (in
  column-major order), but some elements may not be accessed or can be
  restricted to certain values (eg the upper part of a
  lower-triangular matrix)."))

(defclass dense-matrix (dense-matrix-like)
  ()
  (:documentation "Dense matrix, elements stored in column-major order."))

(declaim (inline cm-index2))
(defun cm-index2 (nrow row col)
  "Calculate column-major index, without error checking.  Inlined."
  (+ (* nrow col) row))

;;; xref interface

(defmethod xref ((matrix dense-matrix) &rest subscripts)
  (bind (((row col) subscripts))
    (with-slots (nrow ncol data) matrix
      (check-index row nrow)
      (check-index col ncol)
      (xref data (cm-index2 nrow row col)))))

(defmethod (setf xref) (value (matrix dense-matrix) &rest subscripts)
  (bind (((row col) subscripts))
    (with-slots (nrow ncol data) matrix
      (check-index row nrow)
      (check-index col ncol)
      (setf (xref data (cm-index2 nrow row col))
	    value))))

(defun xsimilar% (object rank)
  ;; not exported
  (let ((lla-type (lla-type object)))
    (ecase rank
      (1 (list (numeric-vector-class lla-type)))
      (2 `(dense-matrix :lla-type ,lla-type)))))

(defmethod xsimilar ((nv numeric-vector) new-dimensions)
  (xsimilar% nv new-dimensions))

(defmethod xsimilar ((m dense-matrix-like) new-dimensions)
  (xsimilar% m new-dimensions))


;;;; matrix creation

(defun initial-contents-from-row-major-list (nrow ncol lla-type elt-list)
  "Return a numeric-vector which has initial contents derived from a flat
list of elements, given in row-major order."
  (check-type nrow dimension)
  (check-type ncol dimension)
  (let* ((size (* nrow ncol))
	 (lisp-type (lla-type->lisp-type lla-type))
	 (data (make-array size :element-type lisp-type)))
    (iter
      (for index :from 0)
      (for elt :in elt-list)
      (bind (((:values row col) (floor index ncol))) ; row major
	(setf (aref data (cm-index2 nrow row col))   ; column major
	      (coerce elt lisp-type))))
    (make-nv size lla-type data t)))

(defgeneric make-matrix (class nrow ncol &key
                               lla-type initial-contents &allow-other-keys)
  (:documentation "Return a matrix of given class and dimensions.
Unless specified otherwise, lla-type is :double and the elements are
initialized with zeros.  Initial contents can be a list or a vector,
in the latter case it may share structure.  Other keys may be allowed
to specify behavior for certain subtypes."))

(defmethod make-matrix (class nrow ncol &key (lla-type :double)
                        (initial-contents 0) (use-directly-p nil))
  "For matrices which are like dense matrices (with perhaps not all of
the elements used, but this is not checked).  For the semantics of
use-directly-p, see make-nv."
  (make-instance class
		 :nrow nrow
		 :ncol ncol
		 :data (etypecase initial-contents
                         (numeric-vector (if (eq (lla-type initial-contents)
                                                 lla-type)
                                             (nv-copy initial-contents)
                                             (nv-copy-convert initial-contents
                                                              lla-type)))
			 ((or null atom) (make-nv (* nrow ncol) lla-type
                                                  initial-contents
                                                  use-directly-p))
			 (list (initial-contents-from-row-major-list
				nrow ncol lla-type initial-contents)))))

(defmethod take ((class (eql 'dense-matrix)) object &key force-copy-p options)
  (declare (ignore force-copy-p))
  (bind (((&key (lla-type :double)) options)
         (dims (xdims object)))
    (unless (= (length dims) 2)
      (error "OBJECT needs to have 2 dimensions for conversion into matrix."))
    (bind (((nrow ncol) dims)
           (contents (make-array (* nrow ncol)
                                 :element-type (lla-type->lisp-type lla-type)))
           (index 0))
      (dotimes (col ncol)
        (dotimes (row nrow)
          (setf (aref contents index) (xref object row col))
          (incf index)))
      (make-matrix nrow ncol lla-type :initial-contents contents
                   :type 'dense-matrix))))

(defmethod xcreate ((class (eql 'dense-matrix)) dimensions &optional options)
  (unless (= (length dimensions) 2)
    (error "Exactly 2 dimensions are needed for a matrix."))
  (bind (((nrow ncol) dimensions)
         ((&key (lla-type :double)) options))
    (make-matrix class nrow ncol :lla-type lla-type)))

(defun matrix-from-first-rows (nv m nrhs n)
  "Extract & return (as a dense-matrix) the first n rows of an m x
nrhs matrix, stored in nv in column-major view.  NOTE: needed to
interface to LAPACK routines like xGELS."
  (let* ((data (nv-data nv))
         (result (make-array (* n nrhs) :element-type (array-element-type data))))
    (dotimes (col nrhs)
      (iter
        (repeat n)
        (for data-index :from (* col m))
        (for result-index :from (* col n))
        (setf (aref result result-index) (aref data data-index))))
    (make-matrix 'dense-matrix n nrhs :lla-type (lla-type nv)
                 :initial-contents result)))

(defun sum-last-rows (nv m nrhs n)
  "Sum & return (as a numeric-vector of the appropriate type) the last
m-n rows of an m x nrhs matrix, stored in nv in column-major view.
NOTE: needed to interface to LAPACK routines like xGELS."
  (let* ((data (nv-data nv))
         (lisp-type (array-element-type data))
         (result (make-array nrhs :element-type lisp-type
                             :initial-element (coerce 0 lisp-type))))
    (dotimes (col nrhs)
      (setf (aref result col)
            (coerce 
             (iter
               (repeat (- m n))
               (for data-index :from (+ (* col m) n))
               (summing (expt (abs (aref data data-index)) 2)))
             lisp-type)))
    (make-nv nrhs (lla-type nv) result t)))

;;;; matrix<->vector conversions
;;;;
;;;; IMPORTANT: results MAY share data with the original.  When
;;;; vectors are converted matrices, they are dense.

(defun vector->matrix-col (numeric-vector)
  "Return vector as a column vector (ie a nx1 dense matrix) of the same
type.  Vector is copied with copy-nv."
  (check-type numeric-vector numeric-vector)
  (make-instance 'dense-matrix :nrow (xsize numeric-vector) :ncol 1
		 :data (nv-copy numeric-vector)))

(defun vector->matrix-row (numeric-vector)
  "Return vector as a row vector (ie a 1xn dense matrix) of the same
type.  Vector is copied with copy-nv."
  (check-type numeric-vector numeric-vector)
  (make-instance 'dense-matrix :nrow 1 :ncol (xsize numeric-vector) 
		 :data (nv-copy numeric-vector)))

(define-condition matrix-not-column-or-row-vector (error)
  ())

(defgeneric matrix->vector (matrix)
  (:documentation "Return an 1xn or nx1 matrix as a numeric-vector of
length n.")
  (:method (matrix)
    ;; fallback method, relies on xref
    (bind (((:slots-read-only nrow ncol) matrix)
	   (vector (make-nv (max nrow ncol) (lla-type matrix))))
      (cond
	((= 1 nrow)			; column matrix
	 (dotimes (i ncol)
	   (setf (xref vector i) (xref 1 i))))
	((= 1 ncol)			; row matrix
	 (dotimes (i ncol)
	   (setf (xref vector i) (xref i 1))))
	(t (error 'matrix-not-column-or-row-vector)))))
  (:method ((matrix dense-matrix))
    ;; fast & cheap method for dense matrices
    (bind (((:slots-read-only nrow ncol data) matrix))
      (unless (or (= 1 nrow) (= 1 ncol))
	(error 'matrix-not-column-or-row-vector))
      (nv-copy data))))

;;;;
;;;;  Semantics of matrices that are subtypes of dense-matrix-like
;;;;  other than dense-matrix: these matrices store their data as if
;;;;  they were dense matrices, but some elements are not accessed and
;;;;  are treated as zero (or some other constant) instead.
;;;;  Nevertheless, these matrices may be treated as dense-matrices,
;;;;  in which case set-restricted can be used to make these elements
;;;;  0.  zero-restricted keeps track of whether this has been done.
;;;;

(define-abstract-class restricted-elements ()
  ((restricted-set-p :type boolean :accessor restricted-set-p
   :initform nil :documentation "non-nil iff not accessible elements
   in the data vectors have been enforced to contain 0."))
  (:documentation "Matrices (& other data structures) that store data
  as some supertype but treat some elements as resticted (eg 0) should
  be subtypes of this, and use restricted-set-p to keep track of
  whether these elements have been set to their values."))

(defgeneric set-restricted* (matrix)
  (:documentation "Always set the restricted elements of matrix,
  regardless of restricted-set-p.  Should not be called directly,
  use set-restricted instead."))

(defgeneric set-restricted (matrix)
  (:documentation "Set restricted/unaccessed elements to the
  appropriate value in the data vector of matrix.  Useful when calling
  functions which expect eg a dense matrix.  If restricted-set-p, or
  if the object is not a subtype of restricted-elements, do nothing.
  NOTE: the default behavior is to call set-restricted* if (not
  restricted-set-p), so it is advised to define that method instead.")
  (:method (matrix))  ;; do nothing
  (:method ((matrix restricted-elements))
    (with-slots (restricted-set-p) matrix
      (unless restricted-set-p
        (set-restricted* matrix)
        (setf restricted-set-p t)))))

(define-abstract-class dense-matrix-restricted (dense-matrix-like restricted-elements)
  ()
  (:documentation "Superclass for restricted dense matrices."))

(defun take-dense-matrix-like (class matrix force-copy-p lla-type)
  "INTERNAL function.  Use data in a dense-matrix matrix in another
  dense-matrix-like object (with given class), copying and converting
  if necessary."
  (with-slots (nrow ncol data) matrix
    (make-instance class :nrow nrow :ncol ncol
                   :data (if (or force-copy-p (not (eq lla-type (lla-type matrix))))
                             (nv-copy-convert data lla-type)
                             (nv-copy data)))))

(defmacro define-dense-matrix-like-take (class)
  "Define a take method for class from dense-matrices."
  `(defmethod take ((class (eql ',class)) (matrix dense-matrix-like)
                    &key force-copy-p options)
     ;; set zeros
     (set-restricted matrix)
     ;; create new structure
     (bind (((&key (lla-type (lla-type matrix))) options))
       (take-dense-matrix-like ',class matrix force-copy-p lla-type))))

(define-dense-matrix-like-take dense-matrix)

;;;;
;;;;  Upper triangular matrix
;;;;
;;;;  (setf xref) has no effect on restricted-set-p.

(defclass upper-triangular-matrix (dense-matrix-like restricted-elements)
  ()
  (:documentation "A dense, upper triangular matrix.  The elements
below the diagonal are not necessarily initialized and not accessed."))

(defmethod set-restricted* ((matrix upper-triangular-matrix))
  (bind (((:slots-read-only nrow ncol data) matrix)
         (vector (copy-data data))
         (zero (coerce 0 (array-element-type vector))))
    ;; set the lower triangle (below diagonal) to 0
    (dotimes (col ncol)
      (iter
        (for index
          :from (1+ (cm-index2 nrow col col))
          :below (cm-index2 nrow nrow col))
        (setf (aref vector index) zero))))
  matrix)

(define-dense-matrix-like-take upper-triangular-matrix)

(defmethod xref ((matrix upper-triangular-matrix) &rest subscripts)
  (bind (((row col) subscripts))
    (with-slots (nrow ncol data) matrix
      (check-index row nrow)
      (check-index col ncol)
      (if (<= row col)
          (xref data (cm-index2 nrow row col))
          (coerce 0 (xelttype matrix))))))

(defmethod (setf xref) (value (matrix upper-triangular-matrix) &rest subscripts)
  (bind (((row col) subscripts))
    (with-slots (nrow ncol data) matrix
      (check-index row nrow)
      (check-index col ncol)
      (if (<= row col)
          (setf (xref data (cm-index2 nrow row col))
                value)
          (unless (zerop value)
            (error 'xref-setting-readonly))))))

(define-matrix-like-printer (matrix upper-triangular-matrix stream) ()
  (print-matrix matrix stream :mask #'upper-triangular-mask))

;;;;
;;;;  Lower triangular matrix
;;;;
;;;;  (setf xref) has no effect on restricted-set-p.

(defclass lower-triangular-matrix (dense-matrix-like restricted-elements)
  ()
  (:documentation "A dense, lower triangular matrix.  The elements
above the diagonal are not necessarily initialized and not accessed."))

(defmethod set-restricted* ((matrix lower-triangular-matrix))
  (bind (((:slots-read-only nrow ncol data) matrix)
         (vector (copy-data data))
         (zero (coerce 0 (array-element-type vector))))
    ;; set the upper triangle (above diagonal) to 0
    (dotimes (col ncol)
      (iter
        (for index
          :from (cm-index2 nrow 0 col)
          :below (cm-index2 nrow col col))
        (setf (aref vector index) zero))))
  matrix)

(define-dense-matrix-like-take lower-triangular-matrix)

(defmethod xref ((matrix lower-triangular-matrix) &rest subscripts)
  (bind (((row col) subscripts))
    (with-slots (nrow ncol data) matrix
      (check-index row nrow)
      (check-index col ncol)
      (if (>= row col)
          (xref data (cm-index2 nrow row col))
          (coerce 0 (xelttype matrix))))))

(defmethod (setf xref) (value (matrix lower-triangular-matrix) &rest subscripts)
  (bind (((row col) subscripts))
    (with-slots (nrow ncol data) matrix
      (check-index row nrow)
      (check-index col ncol)
      (if (>= row col)
          (setf (xref data (cm-index2 nrow row col))
                value)
          (unless (zerop value)
            (error 'xref-setting-readonly))))))

(define-matrix-like-printer (matrix lower-triangular-matrix stream) ()
  (print-matrix matrix stream :mask #'lower-triangular-mask))


;;;;  ??? note on symmetric and hermitian matrices -- Tamas:
;;;;
;;;;  Currently, hermitian matrices can be real, in which case they
;;;;  are also symmetric.  Is this a feature, or lack of good design?
;;;;  Maybe how we handle this could change in the future.  For the
;;;;  present, having a "real hermitian" matrix will just keep this
;;;;  property in operations (eg update-syhe), and functions that
;;;;  accept either a symmetric real or a hermitian complex matrix
;;;;  should only accept a hermitian-matrix, which enforces these
;;;;  properties regardless of whether elements are real or complex.

;;;;
;;;;  Symmetric matrix
;;;;
;;;;  restricted-set-p has to be unset when (setf xref) is used.

(declaim (inline cm-index-symmetric))
(defun cm-index-symmetric (nrow row col)
  (if (<= row col)
      (cm-index2 nrow row col)          ; upper triangle
      (cm-index2 nrow col row)))        ; lower triangle

(defclass symmetric-matrix (dense-matrix-like restricted-elements)
  ()
  (:documentation "A dense symmatric matrix, with elements stored in
   the upper triangle."))

(defmethod xref ((matrix symmetric-matrix) &rest subscripts)
  (bind (((row col) subscripts))
    (with-slots (nrow ncol data) matrix
      (check-index row nrow)
      (check-index col ncol)
      (xref data (cm-index-symmetric nrow row col)))))

(defmethod (setf xref) (value (matrix symmetric-matrix) &rest subscripts)
  (bind (((row col) subscripts))
    (with-slots (nrow ncol data restricted-set-p) matrix
      (check-index row nrow)
      (check-index col ncol)
      (setf restricted-set-p nil)
      (setf (xref data (cm-index-symmetric nrow row col))
	    value))))

(defmethod set-restricted* ((matrix symmetric-matrix))
  (bind (((:slots-read-only nrow ncol data) matrix)
         (vector (copy-data data)))
    ;; set the lower triangle (below diagonal) to mirror elements in
    ;; the upper triangle
    (dotimes (col ncol)
      (iter
        (for row :from col :below nrow)
        (for index
          :from (cm-index2 nrow col col)
          :below (cm-index2 nrow nrow col))
        (setf (aref vector index)
              (aref vector (cm-index2 nrow col row))))))
  matrix)

(define-dense-matrix-like-take symmetric-matrix)

(defun symmetric-mask (row col)
  ;; print = instead of . in the lower triangle
  (if (<= row col)
      nil
      "="))

(define-matrix-like-printer (matrix symmetric-matrix stream) ()
  (print-matrix matrix stream :mask #'symmetric-mask))


;;;;
;;;;  Hermitian matrix
;;;;
;;;;  restricted-set-p has to be unset when (setf xref) is used.

(defclass hermitian-matrix (dense-matrix-like restricted-elements)
  ()
  (:documentation "A dense Hermitian matrix, with elements stored in the upper triangle."))

(defmethod xref ((matrix hermitian-matrix) &rest subscripts)
  (bind (((row col) subscripts))
    (with-slots (nrow ncol data) matrix
      (check-index row nrow)
      (check-index col ncol)
      (if (<= row col)
          (xref data (cm-index2 nrow row col))
          (conjugate (xref data (cm-index2 nrow col row)))))))

(defmethod (setf xref) (value (matrix symmetric-matrix) &rest subscripts)
  (bind (((row col) subscripts))
    (with-slots (nrow ncol data restricted-set-p) matrix
      (check-index row nrow)
      (check-index col ncol)
      (setf restricted-set-p nil)
      (if (<= row col)
          (setf (xref data (cm-index2 nrow row col)) value)
          (setf (xref data (cm-index2 nrow col row)) (conjugate value))))))

(defmethod set-restricted* ((matrix hermitian-matrix))
  (bind (((:slots-read-only nrow ncol data) matrix)
         (vector (copy-data data)))
    ;; set the lower triangle (below diagonal) to conjugate of the
    ;; elements in the upper triangle
    (dotimes (col ncol)
      (iter
        (for row :from col :below nrow)
        (for index
          :from (cm-index2 nrow col col)
          :below (cm-index2 nrow nrow col))
        (setf (aref vector index)
              (conjugate (aref vector (cm-index2 nrow col row)))))))
  matrix)

(define-dense-matrix-like-take hermitian-matrix)

(defun hermitian-mask (row col)
  ;; print * instead of . in the lower triangle
  (if (<= row col)
      nil
      "*"))

(define-matrix-like-printer (matrix hermitian-matrix stream) ()
  (print-matrix matrix stream :mask #'hermitian-mask))


;;;; matrix factorizations
;;;;
;;;; Technically these are not matrices (as they may contain 2+
;;;; matrices and other data), but they are usually defined as
;;;; subclasses of dense-matrix-like (or a subclass, etc) so that you can
;;;; manipulate them as if they were.

(define-abstract-class matrix-factorization ()
  ((nrow :type dimension :initarg :nrow :reader nrow
	 :documentation "The number of rows in the original matrix.")
   (ncol :type dimension :initarg :ncol :reader ncol
	 :documentation "The number of columns in the original matrix.")))

(defgeneric factorization-component (mf component &optional force-copy-p)
  ;;; ??? could lose the force-copy-p with the new nv-copy lazy semantics -- Tamas
  (:documentation "Return a given component of a matrix factorization.
Unless force-copy-p, it can share structure with the original."))

(defgeneric reconstruct (mf)
  (:documentation "Return the original matrix from the matrix factorization."))
  
;;;;
;;;;  LU factorization
;;;;

(defclass LU (matrix-factorization dense-matrix-like)
  ;; LU is not a subclass of matrix, as it is technically two matrices
  ;; (+ pivot indices)
  ;;
  ;; ?? different naming scheme for special matrix types, eg lu would
  ;; become lu-dense, and there would be a superclass named lu? not
  ;; needed at the moment. -- Tamas
  ((ipiv :type numeric-vector-integer :initarg :ipiv :reader ipiv
	 :documentation "pivot indices")
   (data :type numeric-vector :initarg :data
	 :reader data
	 :documentation "LU decomposition of the matrix A with pivoting.")))



;;;; ?? define xarray interface for LU, remember that it has 3
;;;; dimensions, and it is pretty meaningless to access L and U
;;;; separately.  currently, I don't see a need for this, it doesn't
;;;; integrate nicely into xarray.  Let's keep this semi-opauque
;;;; (those who really want it can access ipiv and data and interpret
;;;; it according to LAPACK conventions). -- Tamas

;;;; !! define an initialize-instance after method for consistency check?

;;;; !! define methods for factorization components

;;;;
;;;;  QR factorization
;;;;

(defclass QR (matrix-factorization dense-matrix-like)
  ((data :type numeric-vector :initarg :data
	 :reader data
	 :documentation "QR decomposition of a matrix")))

(defmethod factorization-component ((mf QR) (component (eql :R)) 
                                    &optional force-copy-p)
  (declare (ignore force-copy-p))
  (bind (((:slots-read-only nrow ncol data) mf))
    (take 'upper-triangular-matrix
          (matrix-from-first-rows data nrow ncol ncol))))

;;;; !! define an initialize-instance after method for consistency check?


;;;;
;;;;  Cholesky decomposition
;;;;

(defclass cholesky (matrix-factorization lower-triangular-matrix)
  ((data :type numeric-vector :initarg :data
	 :reader data
	 :documentation "Cholesky decomposition R of a matrix=RR^*,
	 where * is the conjugate transpose.  Should behave as a
	 matrix, but only the lower-triangular part of the
	 decomposition is stored.")))

(defmethod initialize-instance :after ((mf cholesky) &key &allow-other-keys)
  (with-slots (nrow ncol data) mf
    (assert (= nrow ncol))
    (assert (equalp (xdims data) (list (* nrow ncol)))))
  mf)

(defmethod factorization-component ((mf cholesky) component &optional force-copy-p)
  (declare (ignore force-copy-p))
  (take 'lower-triangular-matrix mf))

(define-dense-matrix-like-take cholesky)
