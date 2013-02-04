;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; -*-
(cl:defpackage #:lla-tests
  (:use #:cl
        #:alexandria
        #:anaphora
        #:cl-num-utils
        #:cl-num-utils.matrix-shorthand
        #:cl-slice
        #:clunit
        #:cffi
        #:let-plus
        #:lla)
  (:shadowing-import-from #:alexandria #:mean #:variance #:median)
  (:export #:run))

(in-package :lla-tests)

(defsuite tests ())

(defun run (&optional interactive?)
  "Run all the tests for LLA."
  (run-suite 'tests :use-debugger interactive?))

;; support functions

(defun array= (array1 array2)
  "Test that arrays are equal and have the same element type."
  (and (type= (array-element-type array1)
              (array-element-type array2))
       (equalp array1 array2)))

(defun random-array (type &rest dimensions)
  "Random array for testing."
  (aops:generate* type
                  (if (subtypep type 'complex)
                      (lambda () (coerce (complex (random 100)
                                                  (random 100))
                                         type))
                      (lambda () (coerce (random 100) type)))
                  dimensions))

(defmacro with-foreign-temporary-buffer ((pointer size) &body body)
  "Allocate a buffer for SIZE complex doubles."
  `(with-foreign-pointer (,pointer (* ,size (foreign-type-size :double) 2))
     ,@body))
