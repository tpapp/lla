;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla)

;;;; transpose
;;;
;;; Methods should take care of returning the correct result type (eg
;;; lower-triangular into upper-triangular, etc).  Helper function
;;; transpose-elements% and transpose-matrix% can be used for
;;; implementation.

(defun transpose-elements% (nrow ncol source 
                            &optional 
                            (destination (make-similar-vector source)))
  "Transpose elements from SOURCE into DESTINATION, using NROW and
NCOL, starting at DESTINATION-OFFSET in the latter.  Return DESTINATION."
  (declare (optimize speed (safety 0))
           (fixnum nrow ncol))
  (let ((i 0))
    (with-vector-type-declarations
        (source :other-vectors (destination)
                :simple-test (similar-simple-array? source destination)
                :vector? nil)
      (dotimes (row nrow)
        (dotimes (col ncol)
          (setf (row-major-aref destination i)
                (row-major-aref source (cm-index2 nrow row col)))
          (incf i))))
    destination))

(defun transpose-matrix% (matrix transposed-kind)
  (bind (((:slots-r/o elements nrow ncol) matrix))
    (make-matrix% ncol nrow (transpose-elements% nrow ncol elements) :kind transposed-kind)))

(defmethod transpose ((matrix dense-matrix))
  (transpose-matrix% matrix :dense))

(defmethod transpose ((matrix lower-matrix))
  (transpose-matrix% matrix :upper))

(defmethod transpose ((matrix upper-matrix))
  (transpose-matrix% matrix :lower))

(defmethod transpose ((matrix hermitian-matrix))
  (set-restricted matrix)
  (copy-matrix matrix :copy? t))

;;; conjugate transpose

(defun conjugate-transpose-elements% (nrow ncol elements)
  (declare (optimize speed (safety 0))
           (fixnum nrow ncol))
  (let ((transpose (make-similar-vector elements))
        (i 0))
    (with-vector-type-declarations (elements :other-vectors (transpose))
      (dotimes (row nrow)
        (dotimes (col ncol)
          (setf (aref transpose i)
                (conjugate 
                 (row-major-aref elements (cm-index2 nrow row col))))
          (incf i))))
    transpose))

(defun conjugate-transpose-matrix% (matrix transposed-kind)
  (bind (((:slots-r/o elements nrow ncol) matrix))
    (make-matrix% ncol nrow 
                  (conjugate-transpose-elements% nrow ncol elements)
                  :kind transposed-kind)))

(defgeneric conjugate-transpose (matrix)
  (:documentation "Conjugate transpose of a matrix."))

(defmethod conjugate-transpose ((matrix dense-matrix))
  (conjugate-transpose-matrix% matrix :dense))

(defmethod conjugate-transpose ((matrix lower-matrix))
  (conjugate-transpose-matrix% matrix :upper))

(defmethod conjugate-transpose ((matrix upper-matrix))
  (conjugate-transpose-matrix% matrix :lower))

(defmethod conjugate-transpose ((matrix hermitian-matrix))
  (bind (((:slots-r/o elements nrow) matrix))
    (make-matrix% nrow nrow
                  (with-vector-type-declarations (elements)
                    (map (type-of elements) #'conjugate elements))
                  :kind :hermitian)))
