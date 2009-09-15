(in-package #:lla-asd)

(defpackage #:lla
  (:use :common-lisp
	:cl-utilities
        :iterate
        :bind
	:cffi
	:cl-match
	:xarray)
  (:shadowing-import-from :iterate :collecting :collect)
  
  ;;; !!! export stuff. once things stabilize. -- Tamas
)

(defpackage #:lla-unit-tests
  (:documentation "Unit, validation, and regression testing for lla")
  (:use :common-lisp
	:lift
	:lla)
  (:export run-lla-tests))

;;;; ??? !!! maybe an lla-user package for playing around, once things
;;;; stabilize
