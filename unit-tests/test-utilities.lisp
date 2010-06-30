;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla-unit-tests)

(deftestsuite utilities-tests (lla-unit-tests)
  ())

(addtest (utilities-tests)
  approx=-test
  (let ((*allowed-difference* 9))
    (ensure (approx= 1 9))
    (ensure (not (approx= 1 11)))))

(addtest (utilities-tests)
  ==-test
  (ensure (== (clo 1 :/) (clo 1 :/)))
  (ensure-error (== (clo 1 :/) (clo 2 :/))))

(addtest (utilities-tests)
  make-vector-or-matrix-test
  (let ((*lift-equality-test* #'==))
    (ensure-same (make-vector-with-seq 3 :single) (clo :single 0 1 2))
    (ensure-same (make-matrix-with-seq 2 2 :complex-double)
                 (clo :complex-double
                      0 2 :/
                      1 3))))
