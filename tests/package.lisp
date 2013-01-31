(defpackage #:lla-tests
  (:use #:cl
        #:alexandria
        #:cl-num-utils
        #:cffi
        #:iterate
        #:let-plus
        #:anaphora
        #:lift
        #:lla)
  (:shadowing-import-from #:cl-num-utils #:sum) ; also in ITERATE
  (:shadowing-import-from #:alexandria #:mean #:variance #:median)
  (:export #:run))
