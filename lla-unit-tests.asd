(defsystem :lla-unit-tests
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
    :components
    ((:file "setup")
     (:file "utilities")))
   (:module 
    "tests"
    :pathname #P"tests/"
    :components
    ((:file "test-utilities")
     (:file "test-basics")
     (:file "test-pinned-vector")
     (:file "test-operations")
     (:file "test-linear-algebra"))))
  :depends-on
  (:cl-utilities :iterate :metabang-bind :anaphora 
                 :lla :lift :cl-num-utils))
