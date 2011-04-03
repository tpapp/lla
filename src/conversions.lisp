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
           (diagonal-elements (lla-array n (array-lla-type matrix-elements))))
      (dotimes (i n)
        (setf (aref diagonal-elements i)
              (aref matrix-elements (cm-index2 n i i))))
      (make-diagonal% diagonal-elements))))

;;; conversions to matrix

(defgeneric as-matrix (object &key copy? kind &allow-other-keys)
  (:documentation "Convert object to matrix.  Copying is enforced only
  when COPY?, KIND determines the resulting matrix kind."))

(defmethod as-matrix ((matrix dense-matrix-like) &key copy? (kind :dense))
  (copy-matrix matrix :kind kind :copy? copy?))

(defmethod as-matrix ((array array) &key copy? (kind :dense))
  (declare (ignore copy?))
  (bind (((nrow ncol) (array-dimensions array)))
    (make-matrix% nrow ncol (transpose-elements% ncol nrow array)
                  :kind kind)))

(defmethod as-matrix ((diagonal diagonal) &key (nrow (size diagonal))
                      (ncol (size diagonal)) copy? (kind :dense))
  (declare (ignore copy?))
  (let* ((diagonal-elements (elements diagonal)))
    (assert (<= (size diagonal) (min nrow ncol)) () "Too few rows or columns.")
    (aprog1 (make-matrix nrow ncol (array-lla-type diagonal-elements)
                         :kind kind :initial-element 0)
      (fill-elements-using-diagonal (elements it) diagonal-elements 0 nrow))))

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

(defmethod as-array ((matrix dense-matrix-like) &key copy?)
  (declare (ignore copy?))
  (bind (((:slots-r/o elements nrow ncol) matrix)
         (array (make-array (list nrow ncol) 
                            :element-type (array-element-type
                                           elements))))
    (set-restricted matrix)
    (transpose-elements% nrow ncol elements array)
    array))

(defmethod as-array ((diagonal diagonal) &key (nrow (size diagonal))
                     (ncol (size diagonal)) copy?)
  (declare (ignore copy?))
  (let ((diagonal-elements (elements diagonal)))
    (aprog1 (lla-array (list nrow ncol) (array-lla-type (elements diagonal)) 0)
      (fill-elements-using-diagonal it diagonal-elements 0 ncol))))

;;; row and col

(defun as-row (vector &optional copy?)
  "Convert vector to a row matrix."
  (check-type vector vector)
  (make-matrix% 1 (length vector) (maybe-copy-vector vector copy?)))

(defun as-column (vector &optional copy?)
  "Convert vector to a column matrix."
  (check-type vector vector)
  (make-matrix% (length vector) 1 (maybe-copy-vector vector copy?)))

(defmethod reshape ((matrix dense-matrix) dimensions (order (eql :column-major))
                    &key (kind :dense) (copy? nil copy??))
  (bind ((#(nrow ncol) (reshape-calculate-dimensions dimensions
                                                     (length (elements matrix)))))
    (apply #'copy-matrix matrix :kind kind
           :nrow nrow :ncol ncol
           (when copy??
             `(:copy? ,copy?)))))

(defmethod reshape ((matrix dense-matrix) dimensions (order (eql :row-major))
                    &key copy? (kind :dense))
  (declare (ignore copy?))
  (bind ((elements (elements matrix))
         (size (length elements))
         (#(nrow ncol) (reshape-calculate-dimensions dimensions size))
         (result-elements (make-similar-array elements)))
    (assert (= size (* nrow ncol)))
    (set-restricted matrix)
    (with-indexing* ((vector (nrow matrix) (ncol matrix)) matrix-index
                     matrix-next :column-major? t)
      (with-indexing* ((vector nrow ncol) result-index result-next
                       :column-major? t)
        (loop
          repeat size
          do (setf (row-major-aref result-elements result-index)
                   (row-major-aref elements matrix-index))
          until (matrix-next))))
    (make-matrix% nrow ncol result-elements :kind kind)))
