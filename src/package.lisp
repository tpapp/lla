(defpackage #:lla
  (:use common-lisp iterate bind cffi anaphora alexandria cl-num-utils)
  (:shadowing-import-from cl-num-utils mean variance)
  (:export 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; basics
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   ;; utilities
   
   ;; types

   lla-types lla-to-lisp-type atom-lla-type lla-complex? lla-double?
   real-lla-type complex-lla-type zero* coerce* epsilon* lla-array-element-type
   array-manifest-lla-type make-array* convert-lla-array common-lla-type
   pack packf
   
   ;; printing -- nothing is exported

   ;; special-matrices
   
   wrapped-matrix elements make-matrix mref lower-triangular-matrix
   upper-triangular-matrix hermitian-matrix diagonal-matrix
   
   ;; fortran-atoms -- nothing is exported
   
   ;; printing
   
   *print-lla-precision* *pring-matrix-aligned* *print-matrix-padding*
   
   ;; clo
   
   clo

   ;; pinned-vector -- nothing is exported

   ;; factorizations

   reconstruct lu ipiv qr r square-root left-square-root right-square-root cholesky
   root
   
   ;;    matrix-factorization component reconstruct lu lu-matrix ipiv
   ;;    permutations qr qr-matrix cholesky factor hermitian

   ;; linear-algebra
   
   mm mmm lu solve invert least-squares invert-xx logdet det

   ;;    dot norm1 norm2 normsup outer update-hermitian update-hermitian2
   ;;    hermitian eigen 
   ;;    constrained-least-squares svd tr rank matrix-cond
   
   ))
