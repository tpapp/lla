(defsystem :lla-tests
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
  (:iterate :metabang-bind :anaphora :alexandria :cl-num-utils :cffi :lla :lift))
