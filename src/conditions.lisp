;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla)

;;; errors

(define-condition lapack-error (error) ()
  (:documentation "The LAPACK procedure returned a nonzero info code."))

(define-condition lapack-invalid-argument (lapack-error)
  ((position :initarg :position :type fixnum
             :documentation "Position of the illegal argument"))
  (:documentation "An argument to a LAPACK procedure had an illegal
  value.  Generally, this indicates a bug in LLA and should not
  happen."))

(define-condition lapack-failure (lapack-error)
  ((info :initarg :info :type fixnum
          :documentation "INFO corresponding to error message."))
  (:documentation "Superclass of all LAPACK errors with a positive INFO"))

(define-condition lapack-singular-matrix (lapack-failure) ())

(define-condition lla-incompatible-dimensions (lapack-error) ())

;;; efficiency warnings

(define-condition lla-efficiency-warning (warning)
  ())

(defvar *lla-efficiency-warning-array-type* nil
  "Toggles whether arrays with types not recognized as LLA types raise an
LLA-EFFICIENCY-WARNING-ARRAY-TYPE warning.

Effective only when LLA was loaded & compiled with the appropriate settings in
CL-USER::*LLA-CONFIGURATION*.  See the documentation on configuration in the
README.")

(define-condition lla-efficiency-warning-array-type
    (lla-efficiency-warning)
  ((array :initarg :array))
  (:documentation "See *LLA-EFFICIENCY-WARNING-ARRAY-TYPE*."))

(defmethod print-object ((object lla-efficiency-warning-array-type) stream)
  (format stream "Efficiency warning: had to check a ~A elementwise."
          (slot-value object 'array)))

(defvar *lla-efficiency-warning-array-conversion* nil
  "Toggles whether conversion of array elements to another type when used with
foreign functions raises an LLA-EFFICIENCY-WARNING-ARRAY-CONVERSION warning.

Effective only when LLA was loaded & compiled with the appropriate settings in
CL-USER::*LLA-CONFIGURATION*.  See the documentation on configuration in the
README.")

(define-condition lla-efficiency-warning-array-conversion
    (lla-efficiency-warning)
  ((array :initarg :array :documentation "The array that had to be copied.")
   (type :initarg :type :documentation "Required element type.")))
