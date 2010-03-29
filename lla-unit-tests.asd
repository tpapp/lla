(defsystem :lla-unit-tests
  :description "Unit tests for LLA."
  :version "alpha"
  :author "Tamas K Papp <tkpapp@gmail.com"
  :license "Same as LLA--this is part of the LLA library."
  :serial t
  :components
  ((:module 
    "package-init"
    :pathname #P"unit-tests/"
    :components
    ((:file "package")))
   (:module
    "utilities-and-setup"
    :pathname #P"unit-tests/"
    :components
    ((:file "setup")
     (:file "utilities")))
   (:module 
    "tests"
    :pathname #P"unit-tests/"
    :components
    ((:file "test-utilities")
     (:file "test-basics")
     (:file "test-linear-algebra")
     (:file "test-adjustable"))))
  :depends-on
  (:cl-utilities :iterate :metabang-bind :xarray :anaphora :lla
                 :named-readtables :lift))
