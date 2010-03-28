;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla-unit-tests)

(in-readtable lla:v-syntax)

(deftestsuite utilities-tests (lla-unit-tests)
  ())

(addtest (utilities-tests)
  approx=-test
  (let ((*allowed-difference* 9))
    (ensure (approx= 1 9))
    (ensure (not (approx= 1 11)))))

(addtest (utilities-tests)
  ==-test
  (ensure (== #1v(1) #1v(1)))
  (ensure-error (== #1v(1) #1v(2))))

(addtest (utilities-tests)
  make-nv-or-matrix-test
  (let ((*lift-equality-test* #'==))
    (ensure-same (make-nv-with-seq :single 3) #vs(0 1 2))
    (ensure-same (make-matrix-with-seq :complex-double 2 2)
                 #2vcd(0 2 1 3))))
