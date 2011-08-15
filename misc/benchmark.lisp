;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla)

(defun test-mm (&optional (n 1000000))
  (let* ((a (clo :double 1 2 :/ 3 4))
        (product a))
    (dotimes (i n)
      (setf product (mm a product))
      (when (zerop (mod i 100))
        (setf product a)))))

;;; load-time-value: 12.76
;;; without: 24.79
;;; inline lookup: 22.81 
;;; foreign-funcall: 13.218

(time (test-mm))
