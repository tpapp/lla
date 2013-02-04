(defpackage #:lla
  (:use #:common-lisp
        #:alexandria
        #:anaphora
        #:cffi
        #:cl-num-utils
        #:let-plus
        #:iterate)
  (:shadowing-import-from #:cl-num-utils
                          #:mean #:variance #:median ; also in ALEXANDRIA
                          #:sum            ; also in ITERATE
                          )
  (:export
   ;;package-init - no exports
   ;;basics -- no exports
   ;;types
   #:lla-integer
   #:lla-single
   #:lla-double
   #:lla-complex-single
   #:lla-complex-double
   ;;conditions
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
   #:lla-efficiency-warning-array-conversion
   ;;special-matrices
   #:wrapped-matrix
   #:elements
   #:make-matrix
   #:convert-matrix
   #:mref
   #:lower-triangular-matrix
   #:upper-triangular-matrix
   #:hermitian-matrix
   #:diagonal-matrix
   #:make-lower-triangular-matrix
   #:make-upper-triangular-matrix
   #:make-hermitian-matrix
   #:make-diagonal
   #:dense
   #:lower
   #:hermitian
   #:upper
   #:vec
   #:diag
   #:as-matrix
   ;;printing
   #:*print-lla-precision*
   #:*pring-matrix-aligned*
   #:*print-matrix-padding*
   ;;pinned-array - no exports
   ;;factorizations
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
   #:svd-vt
   ;;fortran-call
   #:with-fp-traps-masked
   ;;linear-algebra (note: some symbols already exported from factorizations)
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
