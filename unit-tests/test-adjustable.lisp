;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-
(in-package :lla-unit-tests)

(deftestsuite adjustable-tests (lla-unit-tests)
  ())

(addtest (adjustable-tests)
  adjustable-numeric-vector
  (let ((v (clo 1 2 3 4))
        (anv (make-anv :double 0 :default-expansion 0))
        (*lift-equality-test* #'x=))
    (add anv v)
    (ensure-same v anv)
    (ensure-same 4 (capacity anv))
    (setf (capacity anv t) 2)
    (ensure-same 6 (capacity anv))
    (setf (size anv t) -3)
    (ensure-same anv (clo 1))
    (shrink anv)
    (ensure-same 1 (capacity anv))
    (add anv (clo 7 8 9))
    (add anv 11d0)
    (ensure-same anv (clo 1 7 8 9 11))))

(addtest (adjustable-tests)
  row-adjustable-matrix
  (let ((m (clo 1 2 :/ 3 4))
        (a (clo 5 6))
        (b (clo 7 8 :/ 9 10))
        ;; ncol=1 is deliberate
        (ram (make-ra-matrix :double 0 1))
        (*lift-equality-test* #'x=))
    (add ram m)
    (ensure-same ram m)
    (add ram b)
    (ensure-same ram (stack-vertically m b) :test #'x=)
    (setf (size ram t) (- (xdim b 1)))
    (ensure-same ram m)
    (ensure-same (capacity ram) (+ (xdim m 1) (xdim b 1)))
    (shrink ram)
    (ensure-same ram m)
    (ensure-same (capacity ram) (xdim m 1))
    (add ram a)
    (ensure-same ram (stack-vertically m a))
    (ensure-same (capacity ram) (1+ (xdim m 1)))
    (add ram (elements a))
    (ensure-same ram (stack-vertically m a a))))

(addtest (adjustable-tests)
  row-adjustable-matrix-one-column
  (let ((ram (make-ra-matrix :double 0 0))
        (*lift-equality-test* #'x=))
    (add ram 5)
    (ensure-same ram (clo 5 :/))
    (add ram (clo 1 :/ 2 3 4))
    (ensure-same ram (clo 5 :/ 1 2 3 4))))

(addtest (adjustable-tests)
  row-adjustable-matrix-submatrix
  (let ((m (clo 1 2 :/ 3 4))
        (ram (make-ra-matrix :double 0 0 :capacity 6)))
    (add ram m)
    (ensure-same (capacity ram) 6 :test #'=)
    (ensure-same (size ram) 2)
    (ensure-same (nrow ram) 2)
    (ensure-same (submatrix ram 1 2 0 nil)
                 (clo 3 4 :/) :test #'==)))