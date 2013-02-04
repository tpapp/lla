(defpackage #:lla-tests
  (:use #:cl
        #:alexandria
        #:cl-num-utils
        #:cl-num-utils.matrix-shorthand
        #:cffi
        #:iterate
        #:let-plus
        #:anaphora
        #:lift
        #:lla)
  (:shadowing-import-from #:cl-num-utils #:sum) ; also in ITERATE
  (:shadowing-import-from #:alexandria #:mean #:variance #:median)
  (:export #:run))
