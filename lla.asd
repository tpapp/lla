;;; Copyright Tamas Papp 2010-2011.
;;;
;;; Distributed under the Boost Software License, Version 1.0.  (See
;;; accompanying file LICENSE_1_0.txt or copy at
;;; http://www.boost.org/LICENSE_1_0.txt)
;;;
;;; This copyright notice pertains to all files in this library.

(defsystem lla
  :description "Lisp Linear Algebra"
  :version "0.1"
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
     (:file "configuration-interface")
     (:file "configuration")
     (:file "libraries")))
   (:module
    "basics"
    :pathname #P"src/"
    :depends-on ("package-init")
    :serial t
    :components
    ((:file "types")
     (:file "conditions")
     (:file "printing")
     (:file "special-matrices")))
   (:module
    "pinned-array"
    :pathname #P"src/"
    :depends-on ("basics")
    :serial t
    :components
    ((:file "foreign-memory")
     (:file "pinned-array")))
   (:module
    "linear-algebra"
    :pathname #P"src/"
    :depends-on ("basics" "pinned-array")
    :serial t
    :components
     ((:file "factorizations")
      (:file "fortran-call")
      (:file "linear-algebra"))))
  :depends-on
  (iterate let-plus cffi anaphora alexandria cl-num-utils))

(defsystem lla-tests
  :description "Unit tests for LLA."
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
  (iterate anaphora alexandria cl-num-utils cl-num-utils-tests cffi lla lift))
