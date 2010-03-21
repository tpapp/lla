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

;;;  adjustable-numeric-vector

(defclass adjustable-numeric-vector (adjustable numeric-vector-like)
  ((size :reader size :initarg :size))
  (:documentation "Adjustable numeric vector.  CAPACITY is the LENGTH
  of ELEMENTS."))

(defmethod xdims ((anv adjustable-numeric-vector))
  (list (size anv)))

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
      (copy-elements-into elements lla-type 0
                          new-elements lla-type 0
                          new-size)
      (setf elements new-elements
            size new-size))
    new-capacity))

(defmethod add ((anv adjustable-numeric-vector) (x number))
  (bind (((:slots default-expansion elements size) anv)
         (new-size (1+ size)))
    (unless (<= new-size size)
      (setf (capacity anv) (max new-size (+ size default-expansion))))
    (setf (aref elements size) x)
    (incf size))
  anv)

(defmethod add ((anv adjustable-numeric-vector) (vector numeric-vector))
  (bind (((:slots lla-type default-expansion elements size) anv)
         ((:slots (vector-lla-type lla-type) (vector-elements elements)) vector)
         (vector-length (length vector-elements))
         (new-size (+ size vector-length)))
    (unless (<= new-size size)
      (setf (capacity anv) (max new-size (+ size default-expansion))))
    (copy-elements-into vector-elements vector-lla-type 0
                        elements lla-type size vector-length)
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
    (copy-elements-into anv-elements anv-lla-type 0
                        (elements result) lla-type 0
                        size)
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


;;;  row-adjustable-matrix

(define-dense-matrix-subclass row-adjustable (adjustable)
  "(* CAPACITY NCOL) ELEMENTS are allocated, but only the first NROW
  are filled.  NCOL will automatically be adjusted by ADD when NROW is
  zero."
  ((capacity :reader capacity :initarg :capacity
             :documentation "Also the leading dimension of the matrix.")))

(defmethod leading-dimension ((ram row-adjustable-matrix))
  (capacity ram))

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
                     elements lla-type capacity
                     new-elements lla-type new-capacity)
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
  (bind (((:slots elements nrow ncol default-expansion capacity lla-type) ram)
         ((:slots-read-only (nrow-matrix nrow) (ncol-matrix ncol)) matrix)
         (new-nrow (+ nrow nrow-matrix)))
    (if (zerop nrow)
        (setf (ncol ram) ncol-matrix)
        (assert (= ncol ncol-matrix) () "Columns don't match."))
    (unless (<= new-nrow capacity)
      (setf (capacity ram) (max new-nrow
                                (+ nrow default-expansion))))
    (set-restricted matrix)
    (copy-columns% nrow-matrix ncol
                   (elements matrix) (lla-type matrix) (leading-dimension matrix)
                   elements lla-type capacity nrow)
    (setf nrow new-nrow))
  ram)

(defmethod add ((ram row-adjustable-matrix) (nv numeric-vector))
  (add ram (vector->row nv)))

(defmethod add ((ram row-adjustable-matrix) (diagonal diagonal))
  (add ram (diagonal->matrix diagonal)))

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
  (bind (((:accessors-r/o lla-type elements nrow ncol leading-dimension) ram)
         ((&key (lla-type lla-type)) options)
         (result (make-matrix lla-type nrow ncol)))
    (copy-columns% nrow ncol
                   elements (lla-type ram) leading-dimension
                   (elements result) lla-type nrow)
    result))

