;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla)

;;; conversion to diagonal

(defgeneric as-diagonal (object &optional copy?)
  (:documentation "Convert object to diagonal, or extract its diagonal
  elements.")
  (:method ((vector vector) &optional copy?)
    (make-diagonal% (if copy?
                        (copy-vector vector)
                        (as-simple-array1 vector))))
  (:method ((matrix dense-matrix-like) &optional copy?)
    (declare (ignore copy?))
    (bind (((:slots-read-only nrow ncol (matrix-elements elements)) matrix)
           (n (min nrow ncol))
           (diagonal-elements (lla-vector n (array-lla-type matrix-elements))))
      (dotimes (i n)
        (setf (aref diagonal-elements i)
              (aref matrix-elements (cm-index2 n i i))))
      (make-diagonal% diagonal-elements))))

;;; conversions to matrix

(defgeneric as-matrix (object &key copy? kind &allow-other-keys)
  (:documentation "Convert object to matrix.  Copying is enforced only
  when COPY?, KIND determines the resulting matrix kind."))

(defmethod as-matrix ((array array) &key copy? (kind :dense))
  (declare (ignore copy?))
  (bind (((nrow ncol) (array-dimensions array)))
    (make-matrix% nrow ncol (transpose-elements% ncol nrow array)
                  :kind kind)))

(defmethod as-matrix ((diagonal diagonal) &key copy? (kind :dense))
  (declare (ignore copy?))
  (let* ((diagonal-elements (elements diagonal))
         (n (length diagonal-elements))
         (matrix-elements (make-similar-vector diagonal-elements
                                               (expt n 2))))
    (fill-matrix-elements-using% diagonal matrix-elements 0 n nil)
    (make-matrix% n n matrix-elements :kind kind)))

(defmethod as-matrix ((vector vector) &key copy? (kind :dense) nrow ncol
                      (row-major? t))
  (bind ((length (length vector)))
    (if nrow
        (check-type nrow dimension)
        (setf nrow (floor length ncol)))
    (if ncol
        (check-type ncol dimension)
        (setf ncol (floor length nrow)))
    (check-matrix-dimensions% nrow ncol length)
    (make-matrix% nrow ncol 
                  (if row-major?
                      (transpose-elements% ncol nrow vector)
                      (maybe-copy-vector vector copy?))
                  :kind kind)))

;; conversions to (Lisp) array

(defgeneric as-array (object &key copy?)
  (:documentation "Convert object to a Lisp array."))

(defmethod as-array ((matrix dense-matrix-like) &key copy?)
  (declare (ignore copy?))
  (bind (((:slots-r/o elements nrow ncol) matrix)
         (array (make-array (list nrow ncol) 
                            :element-type (array-element-type
                                           elements))))
    (set-restricted matrix)
    (transpose-elements% nrow ncol elements array)
    array))

(defmethod as-array ((diagonal diagonal) &key copy?)
  (declare (ignore copy?))
  (let* ((diagonal-elements (elements diagonal))
         (n (length diagonal-elements))
         (array (make-similar-array diagonal-elements (list n n) )))
    (fill-matrix-elements-using% diagonal array 0 n nil)
    array))

;;; row and col

(defun as-row (vector &optional copy?)
  "Convert vector to a row matrix."
  (check-type vector vector)
  (make-matrix% 1 (length vector) (maybe-copy-vector vector copy?)))

(defun as-column (vector &optional copy?)
  "Convert vector to a column matrix."
  (check-type vector vector)
  (make-matrix% (length vector) 1 (maybe-copy-vector vector copy?)))

(defun reshape (object nrow ncol &key (kind :dense) (copy? nil copy??))
  "Reshape vector or DENSE-MATRIX-LIKE object or VECTOR (taken to be
column-major) to a matrix."
  (etypecase object
    (vector (make-matrix% nrow ncol 
                          (maybe-copy-vector object copy?)
                          :kind kind))
    (dense-matrix-like
       (apply #'copy-matrix object :kind kind
              :nrow nrow :ncol ncol
              (when copy??
                `(:copy? ,copy?))))))
