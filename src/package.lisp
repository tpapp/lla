(in-package #:lla-asd)

(defpackage #:lla
  (:use :common-lisp :cl-utilities :iterate :bind :cffi :xarray
        :anaphora :named-readtables)
  (:shadowing-import-from :iterate :collecting :collect)
  (:export 
   ;; utilities -- nothing is exported
   
   ;; types
   
   dimension lla-type lla-complex-p lla-double-p not-within-lla-type
   lisp-type->lla-type lla-type->lisp-type *force-float*

   ;; fortran-atoms -- nothing is exported
   
   ;; numeric-vector
   
   numeric-vector elements nv-array-type numeric-vector-integer
   numeric-vector-single numeric-vector-double
   numeric-vector-complex-single numeric-vector-complex-double
   make-nv* make-nv create-nv ensure-unshared copy-elements
   convert-elements copy-nv convert-nv

   ;; numeric-vector-wrappers -- nothing is exported

   ;; matrix-classes

   matrix-storage nrow ncol matrix-storage-square dense-matrix-like
   matrix-factorization matrix-class matrix-kind

   dense-matrix dense-matrix-integer dense-matrix-single
   dense-matrix-double dense-matrix-complex-single
   dense-matrix-complex-double
   
   restricted-elements

   upper-triangular-matrix upper-triangular-matrix-integer
   upper-triangular-matrix-single upper-triangular-matrix-double
   upper-triangular-matrix-complex-single
   upper-triangular-matrix-complex-double

   lower-triangular-matrix lower-triangular-matrix-integer
   lower-triangular-matrix-single lower-triangular-matrix-double
   lower-triangular-matrix-complex-single
   lower-triangular-matrix-complex-double

   hermitian-matrix hermitian-matrix-integer hermitian-matrix-single
   hermitian-matrix-double hermitian-matrix-complex-single
   hermitian-matrix-complex-double
   
   factorization-component reconstruct

   lu-factorization lu-factorization-integer lu-factorization-single
   lu-factorization-double lu-factorization-complex-single
   lu-factorization-complex-double ipiv

   qr-factorization qr-factorization-integer qr-factorization-single
   qr-factorization-double qr-factorization-complex-single
   qr-factorization-complex-double
   
   cholesky-factorization cholesky-factorization-integer cholesky-factorization-single
   cholesky-factorization-double cholesky-factorization-complex-single
   cholesky-factorization-complex-double

   ;; printing
   
   *pring-matrix-aligned* *print-matrix-padding*
   *print-matrix-precision*
   
   ;; matrix-implementation
   
   cm-index2 make-matrix* make-matrix create-matrix copy-matrix
   convert-matrix vector->matrix vector->column vector->row transpose

   ;; readmacros

   lla-readtable

   ;; fortran-call -- nothing is exported

   ;; linear-algebra

   mmm mm lu solve invert eigen least-squares least-squares-xxinverse
   cholesky 
   
   ))

(defpackage #:lla-unit-tests
  (:documentation "Unit, validation, and regression testing for lla")
  (:use :common-lisp
	:lift
	:lla)
  (:export run-lla-tests))

;;;; ??? !!! maybe an lla-user package for playing around, once things
;;;; stabilize
