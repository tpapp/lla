;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla)

(eval-when (:compile-toplevel :load-toplevel :execute)
  (case (cffi:foreign-type-size :long)
    (8 
       (push :int64 *features*))
    (4 
       (push :int32 *features*))))
