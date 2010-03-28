;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla)

;;; Submatrix - a view on contiguous row & column indexes

(defclass dense-submatrix (dense-matrix-like)
  ((offset :accessor offset :initarg :offset)
   (leading-dimension :accessor leading-dimension :initarg :leading-dimension)))

(defun submatrix (matrix row-start row-end col-start col-end &key copy-p)
  (declare (optimize (debug 3) (speed 0)))
  (check-type matrix dense-matrix-like)
  (flet ((calc-index (index total start?)
           ;; calculate new index, converting negative specifications
           (cond
             ((null index)
              (if start?
                  (error "~A is not a valid start index" index)
                  total))
             ((minusp index)
              (aprog1 (+ total index)
                (assert (<= 0 it) () "~A~A gives negative index" total index)))
             (t 
              (assert (<= index total) () "index ~A above ~A" index total)
              index))))
    (bind (((:accessors-r/o elements lla-type nrow ncol offset leading-dimension)
            matrix)
           (row-start (calc-index row-start nrow t))
           (new-nrow (- (calc-index row-end nrow nil) row-start))
           (col-start (calc-index col-start ncol t))
           (new-ncol (- (calc-index col-end ncol nil) col-start))
           (new-offset (+ offset (* leading-dimension col-start) row-start)))
      (assert (plusp new-nrow) () "Resulting number of rows is not positive")
      (assert (plusp new-ncol) () "Resulting number of columns is not positive")
      (set-restricted matrix)
      (if copy-p
          (let* ((result (make-matrix lla-type new-nrow new-ncol)))
            (copy-columns new-nrow new-ncol
                          elements new-offset leading-dimension lla-type
                          (elements result) 0 new-nrow)
            result)
          (make-instance 'dense-submatrix 
                         :leading-dimension leading-dimension
                         :offset new-offset
                         :ncol new-ncol
                         :nrow new-nrow
                         :elements elements
                         :lla-type lla-type)))))
