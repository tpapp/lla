;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla)

(defun sum-elements (lla-type elements start end)
  "Sum of elements in the given range.  ELEMENTS is guaranteed to
conform to LLA-TYPE.  For internal use."
  (declare (optimize speed (safety 0))
           (fixnum start end))
  (expand-for-lla-types (lla-type :prologue (ecase lla-type))
    `(,lla-type
      (let ((sum (zero* ,lla-type)))
          (declare (cl:type ,(lla::nv-array-type lla-type) elements)
                   (cl:type ,(lla-type->lisp-type lla-type) sum))
        (iter
          (for (the fixnum index) :from start :below end)
          (declare (iterate:declare-variables))
          (incf sum (aref elements index)))
        sum))))

(eval-when (:compile-toplevel :load-toplevel :execute)
  (defun squared-lla-type (lla-type)
    "Return the lla-type of (* x (conjugate x), when x is of LLA-TYPE."
    (ecase lla-type
      (:integer :integer)
      (:single :single)
      (:double :double)
      (:complex-single :single)
      (:complex-double :double))))

(defun sum-squared-elements (lla-type elements start end)
  "Sum of the squares of elements in the given range (for complex
values, (* x (conjugate x)).  ELEMENTS is guaranteed to conform to
LLA-TYPE.  For internal use."
  (declare (optimize speed (safety 0))
           (fixnum start end))
  (lla::expand-for-lla-types (lla-type :prologue (ecase lla-type)
                                       :exclude-integer-p t)
    (bind ((result-type (squared-lla-type lla-type)))
      `(,lla-type
        (let ((sum (zero* ,result-type)))
          (declare (cl:type ,(lla::nv-array-type lla-type) elements)
                   (cl:type ,(lla-type->lisp-type result-type) sum))
          (iter
            (for (the fixnum index) :from start :below end)
            (declare (iterate:declare-variables))
            (incf sum ,(if (lla-complex-p lla-type)
                           ;; here we need the square *modulus*
                           `(let ((x (aref elements index)))
                              (+ (expt (realpart x) 2) (expt (imagpart x) 2)))
                           `(expt (aref elements index) 2))))
          sum)))))

(defun subtract-from-elements (lla-type elements start end value)
  "Subtract value from the elements in the given range.  ELEMENTS is
guaranteed to conform to LLA-TYPE.  For internal use."
  (declare (optimize speed (safety 0))
           (fixnum start end))
  (lla::expand-for-lla-types (lla-type :prologue (ecase lla-type))
    `(,lla-type
      (let ((value (coerce* value ,lla-type)))
        (declare (cl:type ,(lla::nv-array-type lla-type) elements)
                 (,(lla-type->lisp-type lla-type) value))
        (iter
          (for (the fixnum index) :from start :below end)
          (declare (iterate:declare-variables))
          (decf (aref elements index) value))))))

