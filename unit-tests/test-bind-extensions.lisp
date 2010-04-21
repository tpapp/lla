;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla-unit-tests)

(deftestsuite bind-extensions-tests (lla-unit-tests)
  ())

(addtest (bind-extensions-tests)
  lla-vector
  (bind (((:lla-vector a :length length) (make-nv :double 5)))
    (dotimes (i length)
      (setf (a i) (coerce* i :double)))
    (dotimes (i length)
      (incf (a i)))
    (ensure-same a (clo 1 2 3 4 5) :test #'==)))

(addtest (bind-extensions-tests)
  lla-matrix
  (bind (((:lla-matrix a :nrow a-nrow :ncol a-ncol) (make-matrix :double 3 4)))
    (dotimes (row a-nrow)
      (dotimes (col a-ncol)
        (setf (a (a-index row col)) (+ (* 10d0 row) col))))
    (ensure-same a (clo 0 1 2 3 :/
                        10 11 12 13
                        20 21 22 23)
                 :test #'==)))
