(in-package :lla)

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

;;; LLA element type, from the type of data
(defmethod nv-element-type ((matrix matrix))
  (nv-element-type (data matrix)))

;;; Some of the xarray interface, rest defined for specific subclasses

(defmethod xtype ((matrix matrix))
  (xtype (data matrix)))

(defmethod xrank ((matrix matrix))
  2)

(defmethod xdims ((matrix matrix))
  (list (nrow matrix) (ncol matrix)))

(defmethod xdim ((matrix matrix) axis-number)
  (ecase axis-number
    (0 (nrow matrix))
    (1 (ncol matrix))))

;;; printing

(defvar *print-matrix-aligned* t
  "If non-nil, characters will be aligned.")

(defvar *print-matrix-paddig* 1
  "Number of spaces between columns.")

(defun print-matrix (matrix stream &key 
		     (formatter #'standard-numeric-formatter))
  "Format and print the elements of matrix to stream, using
*print-matrix-padding* spaces between columns.  If
*print-matrix-aligned*, columns will be right-aligned.  Prints at most
*print-length* rows and columns, indicating more with a ..."
  ;; ?? maybe column & row labels, not a high priority at the moment
  (bind (((:values nrow row-trunc-p) (print-length-truncate (nrow matrix)))
	 ((:values ncol col-trunc-p) (print-length-truncate (ncol matrix)))
	 (formatted-elements (make-array (list nrow ncol)))
	 (column-widths (make-array ncol :element-type 'fixnum :initial-element 0))
	 (padding (make-array *print-matrix-paddig* :element-type 'character :initial-element #\space))
	 (aligned-p *print-matrix-aligned*))
    ;; first pass - format elements, measure width
    (dotimes (col ncol)
      (dotimes (row nrow)
	(let* ((formatted-element (funcall formatter (xref matrix row col)))
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

(defmethod print-object ((matrix matrix) stream)
  "Uses only the xref for element access.  Respects *print-length* and
aligns output."
    ;; print
    (print-unreadable-object (matrix stream :type t)
      ;; print dimensions and numeric type, as we are not dispatching
      ;; on it here
      (format stream "~a x ~a ~a elements~&" (nrow matrix) (ncol matrix)
	      (nv-element-type matrix))
      ;; print elements
      (print-matrix matrix stream)))

;;; for xdims*, xsize, and xref-writable-p, using defaults.

;;; xref and (setf xref) are defined later.

;;;;
;;;;  Dense matrix class
;;;;

(defclass dense-matrix (matrix)
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

;;;; matrix creation

(defun initial-contents-from-row-major-list (nrow ncol lla-type elt-list)
  "Return a numeric-vector which has initial contents derived from a
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
	(setf (aref data (cm-index2 nrow row col))	     ; column major
	      (coerce elt lisp-type))))
    (make-numeric-vector size lla-type data)))

(defun make-matrix (nrow ncol lla-type &key initial-contents
		    (type 'dense-matrix))
  "Create a matrix of the given dimensions and LLA-type.
Initial-contents can be a scalar (will be coerced and repeated) or a
flat list of elements."
  (unless (eq type 'dense-matrix)
    (error "currently only dense matrices are implemented."))
  (make-instance 'dense-matrix
		 :nrow nrow
		 :ncol ncol
		 :data (etypecase initial-contents
			 ((or null atom) (make-numeric-vector 
					  (* nrow ncol) lla-type
					  initial-contents))
			 (list (initial-contents-from-row-major-list
				nrow ncol lla-type initial-contents)))))

(defmethod take ((class (eql 'dense-matrix)) object &key lla-type)
  (bind ((dims (xdims object)))
    (unless (= (length dims) 2)
      (error "OBJECT needs to have 2 dimensions for conversion into matrix."))
    (bind (((nrow ncol) dims)
           (contents (make-array (* nrow ncol) :element-type (lla-type->lisp-type lla-type)))
           (index 0))
      (dotimes (col ncol)
        (dotimes (row nrow)
          (setf (aref contents index) (xref object row col))
          (incf index)))
      (make-matrix nrow ncol lla-type :initial-contents contents :type 'dense-matrix))))

(defmethod xcreate ((class (eql 'dense-matrix)) dimensions &key (lla-type :double))
  (unless (= (length dimensions) 2)
    (error "Exactly 2 dimensions are needed for a matrix."))
  (bind (((nrow ncol) dimensions))
    (make-matrix nrow ncol lla-type :type 'dense-matrix)))

(defun matrix-from-first-rows (nv m nrhs n)
  "Extract & return (as a dense-matrix) the first n rows of an m x
nrhs matrix, stored in nv in column-major view.  NOTE: needed to
interface to LAPACK routines like xGELS."
  (let* ((data (numeric-vector-data nv))
         (result (make-array (* n nrhs) :element-type (array-element-type data))))
    (dotimes (col nrhs)
      (iter
        (repeat n)
        (for data-index :from (* col m))
        (for result-index :from (* col n))
        (setf (aref result result-index) (aref data data-index))))
    (make-matrix n nrhs (nv-element-type nv) :initial-contents result)))
  
;;;; matrix<->vector conversions
;;;;
;;;; IMPORTANT: results MAY share data with the original.  When
;;;; vectors are converted matrices, they are dense.

(defun vector->matrix-col (numeric-vector)
  "Return vector as a column vector (ie a nx1 dense matrix) of the same
type.  Will share structure with the argument."
  (check-type numeric-vector numeric-vector)
  (make-instance 'dense-matrix :nrow (xsize numeric-vector) :ncol 1
		 :data numeric-vector))

(defun vector->matrix-row (numeric-vector)
  "Return vector as a row vector (ie a 1xn dense matrix) of the same
type.  Will share structure with the argument."
  (check-type numeric-vector numeric-vector)
  (make-instance 'dense-matrix :nrow 1 :ncol (xsize numeric-vector) 
		 :data numeric-vector))

(define-condition matrix-not-column-or-row-vector (error)
  ())

(defgeneric matrix->vector (matrix)
  (:documentation "Return an 1xn or nx1 matrix as a numeric-vector of
length n.  May share structure with the argument.")
  (:method (matrix)
    ;; fallback method, relies on xref
    (bind (((:slots-read-only nrow ncol) matrix)
	   (vector (make-numeric-vector (max nrow ncol)
					(nv-element-type matrix))))
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
      data)))

;;;;
;;;;  LU factorizations
;;;;

(defclass lu ()
  ;; LU is not a subclass of matrix, as it is technically two matrices
  ;; (+ pivot indices).
  ;;
  ;; ?? different naming scheme for special matrix types, eg lu would
  ;; become lu-dense, and there would be a superclass named lu? not
  ;; needed at the moment. -- Tamas
  ((nrow :type dimension :initarg :nrow :reader nrow
	 :documentation "The number of rows in the original matrix.")
   (ncol :type dimension :initarg :ncol :reader ncol
	 :documentation "The number of columns in the original matrix.")
   (ipiv :type numeric-vector-integer :initarg :ipiv :reader ipiv
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

;;;;
;;;;  QR factorizations
;;;;

(defclass QR ()
  ;; QR is not a subclass of matrix, as it is technically two matrices
  ;; (+ pivot indices).
  ((nrow :type dimension :initarg :nrow :reader nrow
	 :documentation "The number of rows in the original matrix.")
   (ncol :type dimension :initarg :ncol :reader ncol
	 :documentation "The number of columns in the original matrix.")
   (data :type numeric-vector :initarg :data
	 :reader data
	 :documentation "QR decomposition of the matrix A with pivoting.")))

;;;; !! define an initialize-instance after method for consistency check?

