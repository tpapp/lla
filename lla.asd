;;; Copyright Tamas Papp 2010.
;;;
;;; Distributed under the Boost Software License, Version 1.0.  (See
;;; accompanying file LICENSE_1_0.txt or copy at
;;; http://www.boost.org/LICENSE_1_0.txt)
;;;
;;; This copyright notice pertains to all files in this library.

(defpackage #:lla-asd
  (:use :cl :asdf))

(in-package #:lla-asd)

(defparameter *fasl-directory*
  (make-pathname :directory '(:relative
			      #+sbcl "fasl-sbcl"
			      #+openmcl "fasl-ccl"
			      #+cmu "fasl-cmucl"
			      #+clisp "fasl-clisp"
			      #-(or sbcl openmcl clisp cmucl) "fasl"
			      )))

(defsystem #:lla
  :description "Lisp Linear Algebra"
  :version "alpha"
  :author "Tamas K Papp <tkpapp@gmail.com>"
  :license "Boost Software License - Version 1.0"
  :serial t
  :components 
  ((:module 
    "package-init"
    :pathname #P"src/"
    :components
    ((:file "package")))
   (:module
    "fortran-interface"
    :pathname #P"src/"
    :serial t
    :components
    ((:file "load-libs")
     (:file "fortran-types")
     (:file "blas-cffi")
     (:file "lapack-cffi")))
   (:module 
    "basics"
    :pathname #P"src/"
    :depends-on ("package-init" "fortran-interface")
    :serial t
    :components
    ((:file "utilities")
     (:file "types")
     (:file "fortran-atoms")
     (:file "copy-elements")
     (:file "matrix-classes")
     (:file "diagonal")
     (:file "printing")
     (:file "bind-extensions")
     (:file "clo")))
   (:module
    "operations"
    :pathname #P"src/"
    :depends-on ("basics")
    :serial t
    :components
    ((:file "elementwise-operations")
     (:file "transpose")
     (:file "conversions")
     (:file "matrix-operations")
     (:file "specialized-utilities")
     (:file "sub")))
   (:module
    "pinned-vector"
    :pathname #P"src/"
    :depends-on ("basics")
    :serial t
    :components
    ((:file "pinned-vector")))
   (:module
    "linear-algebra"
    :pathname #P"src/"
    :depends-on ("basics" "pinned-vector")
    :serial t
    :components
     ((:file "factorizations")
      (:file "lapack-blas-call")
      (:file "linear-algebra"))))
  :depends-on
  (:cl-utilities :iterate :metabang-bind :cffi
                 :anaphora :alexandria :tpapp-utils
                 :cl-num-utils))
