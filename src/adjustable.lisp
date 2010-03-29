;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla)

;;;  Adjustable objects
;;;
;;;  Adjustable objects are flexible along one dimension -- eg vectors
;;;  can be expanded dynamically, and matrices can be augmented with
;;;  new rows or columns.  The implementation is designed such that
;;;  element access is O(1), and the actual adjustment is the costly
;;;  step.  The current extent of the adjustable dimension is referred
;;;  to as SIZE, and the maximum size it can hold without reallocating
;;;  is called CAPACITY.
;;;
;;;  It is recommended that CAPACITY is increased by more than the current
;;;  need if more additions are to follow, this is done automatically
;;;  for ADD, using the DEFAULT-EXPANSION slot.
;;;
;;;  Since the adjustable objects are all subclasses of NUMERIC-VECTOR
;;;  and DENSE-MATRIX-LIKE, they can be mapped to memory and accessed
;;;  directly by foreign calls whenever that is supported by the
;;;  implementation.


;;;  General interface

(define-abstract-class adjustable ()
  ((default-expansion :accessor default-expansion :initarg :default-expansion
     :initform 0 :documentation "Default relative expansion."))
  (:documentation "Object with one dimension adjustable.

Adjustable objects are flexible along one dimension -- eg vectors can
be expanded dynamically, and matrices can be augmented with new rows
or columns.  The implementation is designed such that element access
is O(1), and the actual adjustment is the costly step.  The current
extent of the adjustable dimension is referred to as SIZE, and the
maximum size it can hold without reallocating is called CAPACITY.

It is recommended that CAPACITY is increased by more than the current
need if more additions are to follow, this is done automatically for
ADD, using the DEFAULT-EXPANSION slot.

Since the adjustable objects are all subclasses of NUMERIC-VECTOR and
DENSE-MATRIX-LIKE, they can be mapped to memory and accessed directly
by foreign calls whenever that is supported by the implementation."))

(defgeneric size (adjustable)
  (:documentation "The current SIZE of the adjustable dimension."))

(defgeneric (setf size) (new-size adjustable &optional relative-p)
  (:documentation "Change size to value, or relative to the current
  size if RELATIVE-P.  The new size can only be smaller than the
  current one."))

(defgeneric add (adjustable object)
  (:documentation "Add object an adjustable vector or matrix,
  respectively.  If the object can't take any more elements, its
  CAPACITY is expanded automatically by DEFAULT-EXPANSION."))

(defgeneric capacity (adjustable)
  (:documentation "The current capacity along the adjustable
  dimension."))

(defgeneric (setf capacity) (new-capacity adjustable &optional relative-p)
  (:documentation "Change the capacity to VALUE, relatively to the
  current capacity if RELATIVE-P."))

(defgeneric shrink (adjustable)
  (:documentation "Make capacity equal the current size.")
  (:method ((adjustable adjustable))
    (setf (capacity adjustable) (size adjustable))))

(defun ensure-excess-capacity (adjustable excess-capacity)
  "Adjust capacity so that (<= (+ (SIZE ADJUSTABLE)
EXCESS-CAPACITY) (CAPACITY ADJUSTABLE)).  The inequality may be
strict (see default-expansion).  Return new capacity."
  (bind (((:accessors capacity size default-expansion) adjustable)
         (required-capacity (+ size excess-capacity)))
    (if (<= capacity required-capacity)
        (setf capacity (max required-capacity
                            (+ size default-expansion)))
        capacity)))

;;;  adjustable-numeric-vector

(defclass adjustable-numeric-vector (numeric-vector-like adjustable)
  ((size :reader size :initarg :size :type dimension))
  (:documentation "Adjustable numeric vector.  CAPACITY is the LENGTH
  of ELEMENTS."))

(defmethod xdims ((anv adjustable-numeric-vector))
  (list (size anv)))

(defmethod xsize ((anv adjustable-numeric-vector))
  (size anv))

(defun make-anv (lla-type size &key (initial-element 0) (capacity size)
                 (default-expansion 0))
  "Return an ADJUSTABLE-NUMERIC-VECTOR."
  (assert (<= 0 size capacity))
  (make-instance 'adjustable-numeric-vector 
                 :size size
                 :default-expansion default-expansion
                 :elements (make-nv-elements lla-type capacity initial-element)
                 :lla-type lla-type))

(defmethod print-object ((obj adjustable-numeric-vector) stream)
  (print-unreadable-object (obj stream :type t)
    (with-slots (elements size) obj
      (format stream "~s, capacity ~a with ~a elements: " (lla-type obj) 
              (length elements) size)
      (print-nv-elements elements size stream))))

(defmethod (setf size) (new-size (anv adjustable-numeric-vector) &optional relative-p)
  (bind (((:slots size) anv)
         (new-size (relative-value% size new-size relative-p)))
    (assert (<= 0 new-size (capacity anv)))
    (setf size new-size)))

(defmethod capacity ((anv adjustable-numeric-vector))
  (length (elements anv)))

(defmethod (setf capacity) (new-capacity (anv adjustable-numeric-vector)
                            &optional relative-p)
  (bind (((:slots lla-type size elements) anv)
         (new-capacity (relative-value% (length elements) new-capacity relative-p)))
    (assert (<= 0 new-capacity))
    (let ((new-elements (make-nv-elements lla-type new-capacity))
          (new-size (min size new-capacity)))
      (copy-elements new-size elements 0 lla-type
                     new-elements 0)
      (setf elements new-elements
            size new-size))
    new-capacity))

(defmethod add ((anv adjustable-numeric-vector) (x number))
  (bind (((:slots default-expansion elements size) anv))
    (ensure-excess-capacity anv 1)
    (setf (aref elements size) x)
    (incf size))
  anv)

(defmethod add ((anv adjustable-numeric-vector) (vector numeric-vector))
  (bind (((:slots lla-type default-expansion elements size) anv)
         ((:slots (vector-lla-type lla-type) (vector-elements elements)) vector)
         (vector-length (length vector-elements)))
    (ensure-excess-capacity anv vector-length)
    (copy-elements vector-length vector-elements 0 vector-lla-type
                   elements size)
    (incf size vector-length))
  anv)

(defmethod as* ((class (eql 'adjustable-numeric-vector))
                (nv numeric-vector) copy-p options)
  (declare (ignore copy-p))             ; always copy
  (bind (((&key (lla-type (lla-type nv)) (capacity (xsize nv))) options)
         (result (make-anv lla-type 0 :capacity capacity)))
    (add result nv)))

(defmethod as* ((class (eql 'numeric-vector)) (anv adjustable-numeric-vector)
                copy-p options)
  (declare (ignore copy-p))             ; always copy
  (bind (((:slots-read-only (anv-lla-type lla-type) (anv-elements elements) size) anv)
         ((&key (lla-type anv-lla-type)) options)
         (result (make-nv lla-type size)))
    (copy-elements size anv-elements 0 anv-lla-type
                   (elements result) 0 lla-type)
    result))

(defmethod as* ((class (eql 'vector)) (anv adjustable-numeric-vector)
                copy-p options)
  (declare (ignore copy-p))             ; always copy
  (bind (((:slots-read-only lla-type elements size) anv)
         ((&key (element-type (lla-type->lisp-type lla-type))) options)
         (result (make-array size :element-type element-type)))
    (dotimes (i size)
      (setf (aref result i) (aref elements i)))
    result))

(defmethod vector->column ((anv adjustable-numeric-vector) &key copy-p)
  (declare (ignore copy-p))
  (bind (((:slots-read-only size) anv))
    (make-matrix* (lla-type anv) size 1 (copy-nv-elements anv :length size))))

(defmethod vector->row ((anv adjustable-numeric-vector) &key copy-p)
  (declare (ignore copy-p))
  (bind (((:slots-read-only size) anv))
    (make-matrix* (lla-type anv) 1 size (copy-nv-elements anv :length size))))



;;;  row-adjustable-matrix

(defclass row-adjustable-matrix (adjustable dense-matrix-slots numeric-vector-like)
  ((nrow :type dimension :initarg :nrow :reader nrow
	 :documentation "The number of rows in the matrix.")
   (ncol :type dimension :initarg :ncol :reader ncol
	 :documentation "The number of columns in the matrix.")
   (capacity :reader capacity :initarg :capacity :type dimension
             :documentation "Also the leading dimension of the matrix."))
  (:documentation "(* CAPACITY NCOL) ELEMENTS are allocated, but only
  the first NROW are filled.  NCOL will automatically be adjusted by
  ADD when NROW is zero."))

(defmethod xref ((matrix row-adjustable-matrix) &rest subscripts)
  (bind (((row col) subscripts)
         ((:accessors-r/o elements nrow ncol capacity) matrix))
    (check-index row nrow)
    (check-index col ncol)
    (aref elements (cm-index2 capacity row col))))

(defmethod (setf xref) (value (matrix row-adjustable-matrix) &rest subscripts)
  (bind (((row col) subscripts)
         ((:accessors-r/o elements nrow ncol capacity) matrix))
    (check-index row nrow)
    (check-index col ncol)
    (setf (aref elements (cm-index2 capacity row col))
          value)))

(defun make-ra-matrix (lla-type nrow ncol &key
                       (capacity nrow) (default-expansion 0))
  (assert (<= nrow capacity))
  (make-instance 'row-adjustable-matrix
                 :default-expansion default-expansion
                 :capacity capacity
                 :ncol ncol
                 :nrow nrow
                 :elements (make-nv-elements lla-type (* capacity ncol))
                 :lla-type lla-type))

(defmethod print-object ((matrix row-adjustable-matrix) stream)
  (print-unreadable-object (matrix stream :type t)
    (with-slots (lla-type capacity nrow ncol) matrix
      (format stream "~s, capacity ~a with ~a x ~a elements~&" 
              lla-type capacity nrow ncol)
      (print-matrix matrix stream))))

(defmethod size ((ram row-adjustable-matrix))
  (nrow ram))

(defmethod (setf size) (new-size (ram row-adjustable-matrix) &optional relative-p)
  (bind (((:slots nrow) ram)
         (new-size (relative-value% nrow new-size relative-p)))
    (assert (<= 0 new-size (capacity ram)))
    (setf nrow new-size)))

(defmethod (setf capacity) (new-capacity (ram row-adjustable-matrix)
                            &optional relative-p)
  (bind (((:slots lla-type nrow ncol elements capacity) ram)
         (new-capacity (relative-value% capacity new-capacity relative-p)))
    (assert (<= 0 new-capacity))
    (let ((new-elements (make-nv-elements lla-type (* new-capacity ncol)))
          (new-nrow (min nrow new-capacity)))
      (copy-columns% new-nrow ncol 
                     elements 0 capacity lla-type
                     new-elements 0 new-capacity)
      (setf elements new-elements
            nrow new-nrow
            capacity new-capacity))
    new-capacity))

(defmethod (setf ncol) (new-ncol (ram row-adjustable-matrix))
  ;; ?? Is it resctrictive that we can only change ncol when there are
  ;; no elements in the matrix?  I think it is reasonable, otherwise
  ;; we would need semantics for new elements etc.  The purpose of
  ;; this feature is to facilitate automagically guessing NCOL at the
  ;; first ADD operation.
  (bind (((:slots lla-type ncol nrow elements capacity) ram))
    (assert (zerop nrow) () "NCOL can only be changed when (SIZE=)NROW=0.")
    (setf ncol new-ncol)
    (setf elements (make-nv-elements lla-type (* capacity ncol)))
    ncol))

(defmethod add ((ram row-adjustable-matrix) (matrix dense-matrix-like))
  (bind (((:slots elements nrow ncol capacity lla-type) ram)
         ((:slots-read-only (nrow-matrix nrow) (ncol-matrix ncol)) matrix))
    (if (zerop nrow)
        (setf (ncol ram) ncol-matrix)
        (assert (= ncol ncol-matrix) () "Columns don't match."))
    (set-restricted matrix)
    (ensure-excess-capacity ram nrow-matrix)
    (copy-columns% nrow-matrix ncol
                   (elements matrix) 0 nrow-matrix (lla-type matrix)
                   elements nrow capacity)
    (incf nrow nrow-matrix))
  ram)

(defmethod add ((ram row-adjustable-matrix) (nv numeric-vector))
  (add ram (vector->row nv)))

(defmethod add ((ram row-adjustable-matrix) (diagonal diagonal))
  (add ram (diagonal->matrix diagonal)))

(defmethod add ((ram row-adjustable-matrix) (x number))
  ;; only works for matrices of 1 column
  (bind (((:slots elements nrow ncol lla-type) ram))
    (if (zerop nrow)
        (setf (ncol ram) 1)
        (assert (= ncol 1) () "Columns don't match."))
    (ensure-excess-capacity ram 1)
    (setf (aref elements nrow) (coerce* x lla-type))
    (incf nrow))
  ram)

(defmethod as* ((class (eql 'row-adjustable-matrix)) (matrix dense-matrix-like)
                copy-p options)
  (declare (ignore copy-p))
  (set-restricted matrix)
  (bind (((:slots-read-only lla-type nrow ncol) matrix)
         ((&key (lla-type lla-type) (capacity nrow)) options)
         (result (make-ra-matrix lla-type 0 ncol :capacity capacity)))
    (add result matrix)))

(defmethod as* ((class (eql 'dense-matrix)) (ram row-adjustable-matrix)
                copy-p options)
  (declare (ignore copy-p))
  (bind (((:accessors-r/o lla-type elements nrow ncol capacity) ram)
         ((&key (lla-type lla-type)) options)
         (result (make-matrix lla-type nrow ncol)))
    (copy-columns% nrow ncol
                   elements 0 capacity (lla-type ram)
                   (elements result) 0 nrow)
    result))

(defmethod submatrix ((matrix row-adjustable-matrix) row-start row-end col-start col-end)
  (declare (optimize (debug 3) (speed 0)))
  (bind (((:accessors-r/o elements lla-type nrow ncol capacity) matrix)
         ((:values new-nrow new-ncol offset)
          (submatrix-index-calculations% nrow ncol row-start row-end
                                         col-start col-end capacity)))
    (let ((result (make-matrix lla-type new-nrow new-ncol)))
      (copy-columns new-nrow new-ncol
                    elements offset capacity lla-type
                    (elements result) 0 new-nrow)
      result)))
