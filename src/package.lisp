(defpackage #:lla
  (:use #:common-lisp
        #:alexandria
        #:anaphora
        #:cffi
        #:cl-num-utils
        #:cl-slice
        #:let-plus)
  (:shadowing-import-from #:cl-num-utils ; also in ALEXANDRIA
                          #:mean #:variance #:median)
  ;; no exports from:
  ;;   configuration-interface
  ;;   configuration
  ;;   libraries
  ;;   foreign-memory
  ;;   pinned-array
  (:export                              ; types

   #:lla-integer
   #:lla-single
   #:lla-double
   #:lla-complex-single
   #:lla-complex-double)
  (:export                              ; conditions
   #:lla-internal-error
   #:lla-unhandled-type
   #:lapack-error
   #:lapack-invalid-argument
   #:lapack-failure
   #:lapack-singular-matrix
   #:lla-incompatible-dimensions
   #:lla-efficiency-warning
   #:*lla-efficiency-warning-array-type*
   #:lla-efficiency-warning-array-type
   #:*lla-efficiency-warning-array-conversion*
   #:lla-efficiency-warning-array-conversion)
  (:export                              ; factorizations
   #:lu
   #:ipiv
   #:qr
   #:qr-r
   #:matrix-square-root
   #:xx
   #:left-square-root
   #:right-square-root
   #:cholesky
   #:hermitian-factorization
   #:spectral-factorization
   #:spectral-factorization-w
   #:spectral-factorization-z
   #:svd
   #:svd-u
   #:svd-d
   #:svd-vt)
  (:export                              ; fortran-call
   #:with-fp-traps-masked)
  (:export                              ; linear-algebra
   #:mm
   #:mmm
   #:outer
   #:solve
   #:invert
   #:least-squares
   #:invert-xx
   #:eigenvalues
   #:logdet
   #:det
   #:tr))
