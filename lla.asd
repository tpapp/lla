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
  :depends-on (#:anaphora
               #:alexandria
               #:cffi
               #:cl-num-utils
               #:cl-slice
               #:let-plus)
  :pathname #P"src/"
  :serial t
  :components
  ((:file "package")
   (:file "configuration-interface")
   (:file "configuration")
   (:file "libraries")
   (:file "conditions")
   (:file "types")
   (:file "foreign-memory")
   (:file "pinned-array")
   (:file "factorizations")
   (:file "fortran-call")
   (:file "linear-algebra")
   (:file "blas")))

(defsystem lla-tests
  :description "Unit tests for LLA."
  :author "Tamas K Papp <tkpapp@gmail.com>"
  :license "Same as LLA--this is part of the LLA library."
  :depends-on (#:lla
               #:clunit)
  :pathname #P"tests/"
  :serial t
  :components
  ((:file "setup")
   (:file "pinned-array")
   (:file "linear-algebra")
   (:file "blas")))
