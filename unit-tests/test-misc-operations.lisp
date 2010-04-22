;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla-unit-tests)

(deftestsuite misc-operations-tests (lla-unit-tests)
  ())

(addtest (misc-operations-tests)
  mean
  (ensure-same (mean (clo 1 2 3)) 2d0)
  (ensure-same (mean (clo :integer 1 8 9)) 6))

(addtest (misc-operations-tests)
  sse
  (ensure-same (sse (clo 1 2 3)) 2d0))
