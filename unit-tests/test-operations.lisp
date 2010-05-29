;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla-unit-tests)

(deftestsuite operations-tests (lla-unit-tests)
  ()
  (:equality-test #'==))

;;; elementwise-operations

(addtest (operations-tests)
  emap
  ;; test operations
  (ensure-same (emap #'+ 1 2 (clo :double 3 4))
               (clo :double 6 7))
  (ensure-same (emap #'* (clo :double 3 4)
                     (clo :single 1 2))
               (clo :double 3 8))
  (ensure-same (e- (clo :integer
                        1 2 :/
                        3 4))
               (clo :integer 
                    -1 -2 :/
                    -3 -4))
  (ensure-same (e* (clo :diagonal 2 3)
                   (clo :diagonal 4 5)
                   3d0)
               (clo :double :diagonal
                    24 45))
  (ensure-same (e/ (clo 1 2
                        3 4)
                   (clo 9 1
                        1 9))
               (clo 1/9 2
                    3 4/9))
  (ensure-same (e- #2A((1 2)
                       (3 4))
                   1d0)
               (as-array (clo :double
                              0 1 :/
                              2 3))
               :test #'equalp)
  ;; test errors
  (ensure-error (emap #'+ (clo 1 2 :/ 3 4) (clo :diagonal 1 2 3)))
  (ensure-error (emap #'+ (clo :diagonal 1 2) (clo :diagonal 1 2 3))))

;;; transpose

(addtest (operations-tests)
  transpose
  (let ((m (clo 1 2 3 :/
                4 5 6))
        (mt (clo 1 4 :/
                 2 5
                 3 6)))
    (ensure-same (transpose m) mt)
    (ensure-same (transpose (clo :hermitian 1 2 :/ * 4))
                 (clo :hermitian 1 2 :/ * 4))
    (ensure-same (transpose (copy-matrix m :kind :lower))
                 (copy-matrix mt :kind :upper))
    (ensure-same (transpose (copy-matrix m :kind :upper))
                 (copy-matrix mt :kind :lower))))

(addtest (operations-tests)
  conjugate-transpose
  (let ((m (clo :complex-double
                1 #C(2 9) 3 :/
                4 5 #C(6 -7)))
        (mt (clo :complex-double
                 1 4 :/
                 #C(2 -9) 5
                 3 #C(6 7))))
    (ensure-same (conjugate-transpose m) mt)
    (ensure-same (conjugate-transpose (clo :hermitian 1 #C(2 9) :/ * 4))
                 (clo :hermitian 1 #C(2 -9) :/ * 4))
    (ensure-same (conjugate-transpose (copy-matrix m :kind :lower))
                 (copy-matrix mt :kind :upper))
    (ensure-same (conjugate-transpose (copy-matrix m :kind :upper))
                 (copy-matrix mt :kind :lower))))

;;; conversions

(addtest (operations-tests)
  as-diagonal-tests
  (let ((d (clo :diagonal 5 7)))
    (ensure-same (as-diagonal 
                  (clo 5 0 0 :/
                       0 7 0))
                 d)
    (ensure-same (as-diagonal (clo 5 7)) d)))

(addtest (operations-tests)
  as-matrix-tests
  (ensure-same (as-matrix (clo :diagonal :double 1 2 3))
               (clo :double
                    1 0 0 :/
                    0 2 0
                    0 0 3))
  (ensure-same (as-matrix #2A((1 2) (3 4)))
               (clo 1 2 :/ 3 4))
  (ensure-same (as-matrix (clo :double 1 2 3 4) :nrow 2
                          :kind :hermitian)
               (clo :double :hermitian
                    1 2 :/
                    * 4))
  (ensure-same (as-matrix (clo :double 1 2 3 4) :nrow 2
                          :row-major? nil)
               (clo :double
                    1 3 :/
                    2 4))
  (ensure-error (as-matrix (clo 1 2 3) :nrow 2))
  (ensure-error (as-matrix (clo 1 2 3))))

(addtest (operations-tests)
  as-array-tests
  (let ((*lift-equality-test* #'equalp))
    (ensure-same (as-array (clo :diagonal 1 2 3))
                 #2A((1 0 0)
                     (0 2 0)
                     (0 0 3)))
    (ensure-same (as-array (clo 4 5 :/
                                6 7))
                 #2A((4 5)
                     (6 7)))))

(addtest (operations-tests)
  row-col-reshape-test
  (let ((m (clo 1 4 :/
                2 5
                3 6)))
    (ensure-same (reshape (clo 1 3 5 :/
                               2 4 6)
                          3 2)
                 m)
    (ensure-same (reshape (clo 1 2 3 4 5 6) 3 2)
                 m)
    (ensure-same (as-row (clo 1 2 3))
                 (clo 1 2 3 :/))
    (ensure-same (as-column (clo 1 2 3))
                 (clo 1 :/ 2 3))))

;;; matrix-operations

(addtest (operations-tests)
  matrix-from-first-rows-test
  (ensure-same (lla::matrix-from-first-rows (clo 1 2 3 4 5 6) 2 2 3)
               (clo 1 4 :/
                    2 5)))

(addtest (operations-tests)
  stack-test
  (ensure-same (stack-vertically 
                (clo 1 2 3)
                (clo 4d0 5 6 :/
                     7 8 9)
                (clo :diagonal -1 -2 -3))
               (clo 1 2 3 :/
                    4d0 5 6
                    7 8 9
                    -1 0 0
                    0 -2 0
                    0 0 -3))
  (ensure-same (stack-horizontally
                (clo 1 2)
                (clo 3 4 :/
                     5 6)
                (clo :diagonal 7 8))
               (clo 1 3 4 7 0 :/
                    2 5 6 0 8))
  (ensure-same (stack-vertically
                (clo 1 2)
                (clo 3 4 :/
                     5 6)
                (clo :diagonal 7 8)
                #2A((9 10) (11 12) (13 14)))
               (clo 1 2 :/
                    3 4
                    5 6
                    7 0
                    0 8
                    9 10
                    11 12
                    13 14)))

(addtest (operations-tests)
  eye-test
  (ensure-same (eye 2) (clo 1 0 :/ 0 1))
  (ensure-same (eye 1 :kind :hermitian) (clo :hermitian 1))
  (ensure-same (eye 2 :kind :upper :initial-element 3)
               (clo :upper 3 0 :/ 0 3))
  (ensure-same (eye 3 :kind :diagonal :initial-element 5 :lla-type :integer)
               (clo :integer :diagonal 5 5 5)))

;;; specialized-utilities

(addtest (operations-tests)
  sum-and-sse-test
  (let ((*lift-equality-test* #'=)
        (m (clo :double
                1 2 3 :/
                4 5 6
                7 8 9))
        (d (clo :diagonal :double 1 2 3)))
    (ensure-same (sum (elements m)) 45d0)
    (ensure-same (sum m) 45d0)
    (ensure-same (sum d) 6d0)
    (ensure-same (mean (elements m)) 5d0)
    (ensure-same (mean m) 5d0)
    (ensure-same (mean d) 2d0)
    (ensure-same (sse (elements m)) 60d0)
    (ensure-same (sse m) 60d0)
    (ensure-same (sse d 0) 14d0)
    ;; test below should come last, destructively modifies m too
    ;; set-restricted is not called by copy-matrix, but should be
    ;; called by mean-elements%
    (ensure-same (mean (copy-matrix m :kind :upper :copy? nil))
                 (/ 26d0 9))))

;;; sub

(addtest (operations-tests)
  sub-test
  ;; sub of matrix
  (let ((m (clo 0 1 2 3 4 5 :/
                6 7 8 9 10 11
                12 13 14 15 16 17)))
    (sub m t 2)
    (ensure-same (sub m t 2) (clo 2 8 14))
    (ensure-same (sub m t '(2 . 3)) (clo 2 :/ 8 14))
    (ensure-same (sub m 1 3) 9 :test #'=)
    (ensure-same (sub m #(2 1 0) '(-1 . 0)) (clo 17 :/ 11 5)))
  ;; (setf sub) of matrix
  (let ((m (clo 0 1 2 3 :/
                4 5 6 7
                8 9 10 11)))
    (setf (sub m -2 -1) 14)
    (ensure-same m (clo 0 1 2 3 :/
                        4 5 6 14
                        8 9 10 11))
    (setf (sub m '(0 . 2) 1)
          (clo 17 19))
    (ensure-same m (clo 0 17 2 3 :/
                        4 19 6 14
                        8 9 10 11))
    (setf (sub m #(-1 -2) (cons 1 0))
          (clo 20 21 22 :/
               23 24 25))
    (ensure-same m (clo 0 17 2 3 :/
                        4 23 24 25
                        8 20 21 22))
    (setf (sub m '(0 . 2) #(0 2 3))
          #2A((30 31 32)
              (33 34 35)))
    (ensure-same m (clo 30 17 31 32 :/
                        33 23 34 35
                        8 20 21 22))))
