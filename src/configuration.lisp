;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla)

;;; ********************************************************************
;;; See the documentation on how to configure libraries.  It is highly
;;; unlikely that you need to change anything in this file.
;;; ********************************************************************

(defun default-libraries ()
  "Return a list of libraries.  The source conditions on the platform, relying
TRIVIAL-FEATURES.  This function is only called when the libraries were not
configured by the user, see the documentation on how to do that."
  #+linux '("libblas.so.3gf" "liblapack.so.3gf")
  #+windows '("libblas.dll" "liblapack.dll")
  #+darwin '("libblas.dylib" "liblapack.dylib")
  #-(or linux windows darwin)
  (error "Could not figure out defaults for libraries on your platform. ~
          Please read documentation for instructions on manual ~
          configuration."))

(eval-when (:compile-toplevel :load-toplevel :execute)
  (set-feature 'lla::int64 (query-configuration :int64 nil))
  (set-feature 'lla::cffi-pinning
               (query-configuration :cffi-pinning #+sbcl nil
                                                  #-sbcl t))
  (set-feature 'lla::debug (query-configuration :debug nil))
  (let+ ((efficiency-warnings (query-configuration :efficiency-warnings nil))
         ((&flet enable-efficiency-warning (feature keyword)
            (set-feature feature (find keyword efficiency-warnings)))))
    (enable-efficiency-warning 'lla::efficiency-warning-array-type
                               :array-type)
    (enable-efficiency-warning 'lla::efficiency-warning-array-conversion
                               :array-conversion)))
