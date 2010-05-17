;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla-unit-tests)

(deftestsuite sub-tests (lla-unit-tests)
  ())

(addtest (sub-tests)
  sub-nv
  (let ((*lift-equality-test* #'==)
        (a (clo 0 1 2 3 4 5)))
    (ensure-same (sub a t) a)
    (ensure-same (sub a (cons 2 -2)) (clo 2 3))
    (ensure-same (sub a #(1 2 -1 -2)) (clo 1 2 5 4))))

(addtest (sub-tests)
  sub-matrix
  (let ((*lift-equality-test* #'==)
        (a (clo 0 1 2 :/
                3 4 5
                6 7 8)))
    (ensure-same (sub a t t) a)
    (ensure-same (sub a '(1 . 0) '(0 . -1)) (clo 3 4 :/
                                                 6 7))
    (ensure-same (sub a t 1) (clo 1 4 7))
    (sub a 1 -1)
    (ensure-same (sub a 1 -1) 5d0 :test #'equal)))
