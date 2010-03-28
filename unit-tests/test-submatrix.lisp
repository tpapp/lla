;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla-unit-tests)

(in-readtable lla:v-syntax)

(addtest (linear-algebra-tests)
  submatrix
  (let* ((matrix #4v(1 2 3 4
                     5 6 7 8
                     9 10 11 12
                     13 14 15 16))
         (submatrix #2v(6 7 10 11))
         (submatrix1 (submatrix matrix 1 3 1 3))
         (*lift-equality-test* #'x=))
    (ensure-same submatrix submatrix1)
    (ensure-same (submatrix matrix 1 -1 1 -1) submatrix1)
    (ensure-same (solve submatrix1 submatrix) #2v(1 0 0 1)
                 :test (x~= 1d-8))))

(addtest (linear-algebra-tests)
  submatrix-copied
  (let* ((matrix #4v(1 2 3 4
                     5 6 7 8
                     9 10 11 12
                     13 14 15 16))
         (submatrix #2v(6 7 10 11))
         (submatrix1 (submatrix matrix 1 3 1 3 :copy-p t))
         (*lift-equality-test* #'x=))
    (ensure-same submatrix submatrix1)
    (ensure-same (submatrix matrix 1 -1 1 -1 :copy-p t) submatrix1)
    (ensure-same (solve submatrix1 submatrix) #2v(1 0 0 1)
                 :test (x~= 1d-8))))
