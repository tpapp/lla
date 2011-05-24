;;; Copyright Tamas Papp 2010-2011.
;;;
;;; Distributed under the Boost Software License, Version 1.0.  (See
;;; accompanying file LICENSE_1_0.txt or copy at
;;; http://www.boost.org/LICENSE_1_0.txt)
;;;
;;; This copyright notice pertains to all files in this library.

(defsystem lla
  :description "Lisp Linear Algebra"
  :version "beta"
  :author "Tamas K Papp <tkpapp@gmail.com>"
  :license "Boost Software License - Version 1.0"
  :serial t
  :components 
  ((:module 
    "package-init"
    :pathname #P"src/"
    :serial t
    :components
    ((:file "package")
     (:file "libraries")))
   (:module
    "foreign-interface"
    :pathname #P"swig/"
    :components
    ((:file "lla-foreign-interface")))
   (:module 
    "basics"
    :pathname #P"src/"
    :depends-on ("package-init" "foreign-interface")
    :serial t
    :components
    ((:file "utilities")
     (:file "types")
     (:file "foreign-atoms")
     (:file "printing")
     (:file "special-matrices")
     (:file "clo")))
   (:module
    "pinned-array"
    :pathname #P"src/"
    :depends-on ("basics")
    :serial t
    :components
    ((:file "pinned-array")))
   (:module
    "linear-algebra"
    :pathname #P"src/"
    :depends-on ("basics" "pinned-array")
    :serial t
    :components
     ((:file "factorizations")
      (:file "lapack-blas-call")
      (:file "linear-algebra"))))
  :depends-on
  (iterate metabang-bind let-plus cffi anaphora alexandria cl-num-utils))

(defsystem lla-tests
  :description "Unit tests for LLA."
  :version "alpha"
  :author "Tamas K Papp <tkpapp@gmail.com"
  :license "Same as LLA--this is part of the LLA library."
  :serial t
  :components
  ((:module 
    "package-init"
    :pathname #P"tests/"
    :components
    ((:file "package")))
   (:module
    "utilities-and-setup"
    :pathname #P"tests/"
    :serial t
    :components
    ((:file "setup")
     (:file "utilities")))
   (:module 
    "tests"
    :pathname #P"tests/"
    :components
    ((:file "basics")
     (:file "pinned-array")
     (:file "linear-algebra"))))
  :depends-on
  (iterate metabang-bind anaphora alexandria cl-num-utils cffi lla lift))
