(defpackage #:lla-tests
  (:use #:cl #:iterate #:let-plus #:anaphora #:cl-num-utils #:cffi #:lift
        #:lla)
  (:import-from #:cl-num-utils-tests #:array=)
  (:shadowing-import-from #:cl-num-utils #:sum) ; also in ITERATE
  (:export #:run))
