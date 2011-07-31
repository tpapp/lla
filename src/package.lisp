(defpackage #:lla
  (:use common-lisp iterate let-plus cffi anaphora alexandria cl-num-utils)
  (:shadowing-import-from cl-num-utils mean variance)
  (:export 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; basics
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   ;; utilities

   ensure-vector ensure-matrix
   
   ;; types

   lla-types lla-to-lisp-type atom-lla-type lla-complex? lla-double?
   real-lla-type complex-lla-type zero* coerce* epsilon* lla-array-element-type
   array-manifest-lla-type make-array* convert-lla-array common-lla-type
   pack packf
   
   ;; special-matrices
   
   wrapped-matrix elements make-matrix convert-matrix mref
   lower-triangular-matrix upper-triangular-matrix hermitian-matrix
   diagonal-matrix make-lower-triangular-matrix make-upper-triangular-matrix
   make-hermitian-matrix make-diagonal
   
   ;; printing
   
   *print-lla-precision* *pring-matrix-aligned* *print-matrix-padding*
   
   ;; clo
   
   clo

   ;; factorizations

   as-matrix lu ipiv qr qr-r matrix-square-root xx left-square-root
   right-square-root cholesky hermitian-factorization 
   
   spectral-factorization spectral-factorization-w spectral-factorization-z
   
   svd svd-u svd-d svd-vt
   
   ;; fortran-call
   
   with-fp-traps-masked
   
   ;; linear-algebra
   
   mm mmm lu solve invert least-squares invert-xx logdet det tr

   ;;    dot norm1 norm2 normsup outer update-hermitian update-hermitian2
   ;;    hermitian eigen 
   ;;    constrained-least-squares svd tr rank matrix-cond
   
   ))
