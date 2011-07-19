;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla)

;;; ********************************************************************
;;; See the documentation on how to configure libraries.  It is highly
;;; unlikely that you need to change anything in this file.
;;; ********************************************************************

(defun query-configuration (indicator &optional default)
  "Return the property for INDICATOR, with an optional default, which can be a
function (called on demand)."
  (let* ((symbol 'cl-user::*lla-configuration*)
         (plist (when (boundp symbol)
                  (symbol-value symbol)))
         (unique (gensym))
         (configuration (getf plist indicator unique)))
    (if (eq configuration unique)       ; only call when necessary
        (if (functionp default)
            (funcall default)
            default)
        configuration)))

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

