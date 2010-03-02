(defpackage #:lla-unit-tests
    (:use #:cl #:cl-utilities #:iterate #:metabang-bind #:cffi #:xarray
          #:anaphora #:named-readtables #:lift #:lla)
  (:shadowing-import-from :iterate :collecting :collect)
  (:import-from lla coercible-pairs-list mem-aref* 
                with-fortran-atoms make-nv-elements
                with-nv-input with-vector-output copy-matrix%)
  (:export run-lla-tests))
