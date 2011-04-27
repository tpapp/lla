;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; -*-

(in-package #:lla)

;;; !! write a bit of documentation

#+lla-use-mkl
(progn
  (define-foreign-library :gomp
    (:unix (:or "libgomp.so.1")))
  (define-foreign-library :iomp
    (:unix (:or "libiomp5.so")))
  (define-foreign-library :mkl-rt
    (t (:default "libmkl_rt")))
  (define-foreign-library :pthread
    (:unix (:or "libpthread.so.0"))
    (t (:default "libpthread")))
  (use-foreign-library :gomp)
  (use-foreign-library :iomp)
  (use-foreign-library :mkl-rt)
  (use-foreign-library :pthread))

#-lla-use-mkl
(progn
  (error "this part is not written yet, patches are welcome"))
