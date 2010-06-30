;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla)

(defun copy-elements (source source-offset destination destination-offset length)
  "Copy LENGTH elements from SOURCE to DESTINATION (both are arrays),
starting at OFFSET (elements are copied using row-major indexing),
coercing elements if necessary."
  (macrolet ((copy-loop (source-element-form)
               `(iter
                  (for source-index :from source-offset
                       :below (+ source-offset length))
                  (for destination-index :from destination-offset)
                  (setf (row-major-aref destination destination-index)
                        ,source-element-form))))
    (let ((destination-type (array-element-type destination)))
      (if (equal destination-type (array-element-type source))
          (with-vector-type-declarations (source 
                                          :other-vectors (destination)
                                          :simple-test (simple-array? destination))
            (copy-loop (row-major-aref source source-index)))
          (muffle-optimization-notes
            (copy-loop (coerce (row-major-aref source source-index)
                               destination-type))))))    
  (values))

(defun copy-vector (vector &optional (lla-type (array-lla-type vector)))
  "Copy a vector, optionally converting to LLA-TYPE."
  (let* ((length (length vector))
         (destination (lla-vector length lla-type)))
    (copy-elements vector 0 destination 0 length)
    destination))

(defun maybe-copy-vector (vector copy? &optional lla-type)
  "Shallow copying unless LLA-TYPE or COPY?.  The result is always a
SIMPLE-ARRAY1."
  (cond
    (lla-type (copy-vector vector lla-type))
    (copy? (copy-vector vector))
    (t (as-simple-array1 vector))))

(declaim (inline cm-index2)
         (ftype (function (fixnum fixnum fixnum) fixnum)))
(defun cm-index2 (nrow row col)
  "Calculate column-major index, without error checking.  Inlined."
  (the fixnum (+ (the fixnum (* nrow col)) row)))

(defun copy-columns (nrow ncol source source-offset source-ld
                      destination destination-offset destination-ld)
  "Copy elements from one matrix to another, respecting leading
dimensions and offsets.  If types are different, elements are coerced.
Return no values.  Accepts all kinds if simple arrays, uses
row-major-aref."
  (declare (optimize (speed 3) (safety 0))
           (fixnum nrow ncol source-offset source-ld
                   destination-offset destination-ld))
  (macrolet ((copy-loop (source-element-form)
               `(iter
                  (for (the fixnum col) :from 0 :below ncol)
                  (declare (iterate:declare-variables))
                  (iter
                    (for row :from 0 :below nrow)
                    (for source-index
                         :from (+ source-offset (cm-index2 source-ld 0 col)))
                    (for destination-index
                         :from (+ destination-offset (cm-index2 destination-ld 0 col)))
                    (declare (iterate:declare-variables)
                             (fixnum row source-index destination-index))
                    (setf (row-major-aref destination destination-index)
                          ,source-element-form)))))
    (let ((destination-type (array-element-type destination)))
      (if (equal destination-type (array-element-type source))
          (with-vector-type-declarations (source 
                                          :other-vectors (destination)
                                          :simple-test (simple-array? destination))
            (copy-loop (row-major-aref source source-index)))
          (muffle-optimization-notes
            (copy-loop (coerce (row-major-aref source source-index)
                               destination-type))))))
  (values))
