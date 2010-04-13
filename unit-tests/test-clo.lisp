;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla-unit-tests)

(deftestsuite clo-tests (lla-unit-tests)
  ())

(addtest (clo-tests)
  clo-vector
  (let ((*lift-equality-test* #'==))
    (ensure-same (clo 1 2 3)
                 (create-nv '(1 2 3) :double))
    (ensure-same (clo :single 1 2 3)
                 (create-nv '(1 2 3) :single))
    (ensure-same (clo :no-coerce 1d0 2d0 3d0)
                 (create-nv '(1 2 3) :double))
    (ensure-error (clo :no-coerce 1 2 3))))

(addtest (clo-tests)
  clo-diagonal
  (let ((*lift-equality-test* #'==))
    (ensure-same (clo :diagonal 1 2 3)
                 (create-diagonal '(1 2 3) :double))
    (ensure-same (clo :single :diagonal 1 2 3)
                 (create-diagonal '(1 2 3) :single))
    (ensure-same (clo :diagonal :no-coerce 1d0 2d0 3d0)
                 (create-diagonal '(1 2 3) :double))
    (ensure-error (clo :diagonal :no-coerce 1 2 3))))

(addtest (clo-tests)
  clo-matrix
  (let ((*lift-equality-test* #'==))
    (ensure-same (clo :dense 1 2 :/ 3 4)
                 (create-matrix 2 '(1 2 3 4)))
    (ensure-same (clo :dense 1 :/ 2 3 4)
                 (create-matrix 1 '(1 2 3 4)))
    (ensure-same (clo :hermitian 1 2 :/ 3 4)
                 (create-matrix 2 '(1 2 3 4) :kind :hermitian))
    (ensure-same (clo :lower-triangular 1 2 :/ 3 4)
                 (create-matrix 2 '(1 2 3 4) :kind :lower-triangular))
    (ensure-same (clo :upper-triangular 1 2 :/ 3 4)
                 (create-matrix 2 '(1 2 3 4) :kind :upper-triangular))))
