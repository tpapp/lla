;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla)

;;;  Adjustable objects
;;;
;;;  Adjustable objects are flexible along one dimension -- eg vectors
;;;  can be expanded dynamically, and matrices can be augmented with
;;;  new rows or columns.  The implementation is designed such that
;;;  element access is O(1), and the actual adjustment is the costly
;;;  step.  Objects have a SIZE, which tells how much the flexible
;;;  dimension can be augmented without rellocating and copying.  It
;;;  is recommended that SIZE is increased by more than the current
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
  (:documentation "LLA objects, with one dimension adjustable."))

(defgeneric size (adjustable)
  (:documentation ""))

(defgeneric add (adjustable object)
  (:documentation "Add object an adjustable vector or matrix,
  respectively.  If the object can't take any more elements, it is
  expanded automatically by DEFAULT-EXPANSION."))

(defgeneric adjust (adjustable new-size &optional relative-p)
  (:documentation "Set the adjustable dimension to NEW-SIZE.  If
  RELATIVE-P, it is added to current size, in this case, negative
  NEW-SIZE is accepted.  If the resulting size cannot accommodate all
  elements, those at the end are discarded.  Note: implementations of
  SHRINK and CLEAR should just calculate the new size and call this
  method."))

(defgeneric shrink (adjustable)
  (:documentation "Make SIZE equal the current value of the adjustable dimension."))

(defgeneric clear (adjustable &optional count)
  (:documentation "Clear the last COUNT elements (default is all), but
  do not adjust size."))

;;;  adjustable-numeric-vector

(defclass adjustable-numeric-vector (adjustable numeric-vector-like)
  ((len :accessor len :initarg :len))
  (:documentation "Adjustable numeric vector.  SIZE is the LENGTH of
  ELEMENTS, and the actual number of elements is LENGTH%"))

(defmethod size ((anv adjustable-numeric-vector))
  (length (elements anv)))

(defmethod xsize ((anv adjustable-numeric-vector))
  (len anv))

(defun make-anv (lla-type length &key (initial-element 0) (size length)
                 (default-expansion 0))
  "Return an ADJUSTABLE-NUMERIC-VECTOR."
  (assert (<= length size) () "SIZE has to be at least as large as LENGTH.")
  (make-instance 'adjustable-numeric-vector 
                 :len length
                 :default-expansion default-expansion
                 :elements (make-nv-elements lla-type size initial-element)
                 :lla-type lla-type))

(defmethod print-object ((obj adjustable-numeric-vector) stream)
  (print-unreadable-object (obj stream :type t)
    (with-slots (elements len) obj
      (let* ((size (length elements)))
	(format stream "~s, size ~a with ~a elements: " (lla-type obj) size len)
        (print-nv-elements elements len stream)))))

(defmethod adjust ((nv adjustable-numeric-vector) new-size &optional relative-p)
  (bind (((:accessors-read-only lla-type size) nv)
         ((:slots elements len) nv)
         (new-size (+ (if relative-p size 0) new-size)))
    (assert (<= 0 new-size) () "Size has to be nonnegative.")
    (bind ((new-elements (make-nv-elements lla-type new-size))
           (new-len (min len new-size)))
      (copy-elements-into elements lla-type 0
                          new-elements lla-type 0
                          new-len)
      (setf elements new-elements
            len new-len)))
  nv)

(defmethod add ((adjustable adjustable-numeric-vector) (x number))
  (bind (((:accessors-read-only size default-expansion) adjustable)
         ((:slots elements len) adjustable)
         (required-size (1+ len)))
    (unless (<= required-size size)
      (adjust adjustable (max required-size (+ size default-expansion))))
    (setf (aref elements len) x)
    (incf len))
  adjustable)

(defmethod add ((adjustable adjustable-numeric-vector) (nv numeric-vector))
  (bind (((:accessors-read-only lla-type size default-expansion) adjustable)
         ((:slots elements len) adjustable)
         ((:slots (nv-lla-type lla-type) (nv-elements elements)) nv)
         (nv-length (length nv-elements))
         (required-size (+ nv-length len)))
    (unless (<= required-size size)
      (adjust adjustable (max required-size (+ size default-expansion))))
    (copy-elements-into nv-elements nv-lla-type 0
                        elements lla-type len nv-length)
    (incf len nv-length))
  adjustable)

(defmethod shrink ((adjustable adjustable-numeric-vector))
  (adjust adjustable (len adjustable)))

(defmethod clear ((adjustable adjustable-numeric-vector) &optional count)
  (with-slots (len) adjustable
    (decf len (aif count
                   it
                   len)))
  adjustable)

(defmethod as* ((class (eql 'adjustable-numeric-vector)) (nv numeric-vector) copy-p options)
  (declare (ignore copy-p))
  (bind (((&key (lla-type (lla-type nv)) (size (xsize nv))) options)
         (result (make-anv lla-type 0 :size size)))
    (add result nv)))

(defmethod as* ((class (eql 'numeric-vector)) (anv adjustable-numeric-vector) copy-p options)
  (declare (ignore copy-p))
  (bind (((:slots-read-only (anv-lla-type lla-type) (anv-elements elements) len) anv)
         ((&key (lla-type anv-lla-type)) options)
         (result (make-nv lla-type len)))
    (copy-elements-into anv-elements anv-lla-type 0
                        (elements result) lla-type 0
                        len)
    result))

;;;  row-adjustable-matrix

(define-dense-matrix-subclass row-adjustable (adjustable)
  "(* SIZE NCOL) ELEMENTS are allocated, but only the
  first NROW are filled."
  ((leading-dimension :accessor leading-dimension :initarg :leading-dimension)))

(defmethod size ((matrix row-adjustable-matrix))
  (leading-dimension matrix))

(defun make-ra-matrix (lla-type nrow ncol &key
                       (size nrow) (default-expansion 0))
  (assert (<= nrow size) () "SIZE has to be at least as large as NROW.")
  (make-instance 'row-adjustable-matrix :default-expansion default-expansion
                 :leading-dimension size
                 :ncol ncol
                 :nrow nrow
                 :elements (make-nv-elements lla-type (* size ncol))
                 :lla-type lla-type))

(defmethod print-object ((matrix row-adjustable-matrix) stream)
  (print-unreadable-object (matrix stream :type t)
    (with-slots (lla-type leading-dimension nrow ncol) matrix
      (format stream "~s, size ~a with ~a x ~a elements~&" 
              lla-type leading-dimension nrow ncol)
      (print-matrix matrix stream))))

(defmethod adjust ((matrix row-adjustable-matrix) new-size &optional relative-p)
  (bind (((:slots-read-only lla-type ncol) matrix)
         ((:slots elements nrow leading-dimension) matrix)
         (new-leading-dimension (+ (if relative-p leading-dimension 0) new-size)))
    (assert (<= 0 new-leading-dimension) () "Leading dimension has to be nonnegative.")
    (bind ((new-elements (make-nv-elements lla-type (* new-leading-dimension ncol)))
           (new-nrow (min nrow new-leading-dimension)))
      (copy-columns% new-nrow ncol 
                     elements lla-type leading-dimension
                     new-elements lla-type new-leading-dimension)
      (setf elements new-elements
            nrow new-nrow
            leading-dimension new-leading-dimension)))
  matrix)

(defmethod add ((adjustable row-adjustable-matrix) (matrix dense-matrix-like))
  (bind (((:slots lla-type elements nrow ncol leading-dimension
                  default-expansion) adjustable)
         ((:slots-read-only (nrow-matrix nrow) (ncol-matrix ncol)) matrix)
         (required-leading-dimension (+ nrow nrow-matrix)))
    (assert (= ncol ncol-matrix) () "Columns don't match.")
    (unless (<= required-leading-dimension leading-dimension)
      (adjust adjustable (max required-leading-dimension 
                              (+ leading-dimension default-expansion))))
    (set-restricted matrix)
    (copy-columns% nrow-matrix ncol
                   (elements matrix) (lla-type matrix) (leading-dimension matrix)
                   elements lla-type leading-dimension nrow)
    (incf nrow nrow-matrix))
  adjustable)

(defmethod add ((adjustable row-adjustable-matrix) (nv numeric-vector))
  (add adjustable (vector->row nv)))

(defmethod add ((adjustable row-adjustable-matrix) (diagonal diagonal))
  (add adjustable (diagonal->matrix diagonal)))

(defmethod shrink ((adjustable row-adjustable-matrix))
  (adjust adjustable (nrow adjustable)))

(defmethod clear ((adjustable row-adjustable-matrix) &optional count)
  (with-slots (nrow) adjustable
    (decf nrow (aif count
                    it
                    nrow)))
  adjustable)

(defmethod as* ((class (eql 'adjustable-row-matrix)) (matrix dense-matrix-like)
                copy-p options)
  (declare (ignore copy-p))
  (set-restricted matrix)
  (bind (((:slots-read-only lla-type nrow ncol) matrix)
         ((&key (lla-type lla-type) (size nrow)) options)
         (result (make-ra-matrix lla-type 0 ncol :size size)))
    (add result matrix)))

(defmethod as* ((class (eql 'dense-matrix)) (matrix row-adjustable-matrix)
                copy-p options)
  (declare (ignore copy-p))
  (bind (((:slots-read-only lla-type elements nrow ncol leading-dimension) matrix)
         ((&key (lla-type lla-type)) options)
         (result (make-matrix lla-type nrow ncol)))
    (copy-columns% nrow ncol
                   elements (lla-type matrix) leading-dimension
                   (elements result) lla-type nrow)
    result))
