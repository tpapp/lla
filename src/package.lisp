(in-package #:lla-asd)

(defpackage #:lla
    (:use :common-lisp :cl-utilities :iterate :bind :cffi :xarray
          :anaphora :named-readtables :tpapp-utils)
  (:shadowing-import-from :iterate :collecting :collect)
  (:export 
   ;; utilities -- nothing is exported
   
   ;; types
   
   dimension lla-type lla-complex-p lla-double-p not-within-lla-type
   lisp-type->lla-type lla-type->lisp-type zero* coerce* epsilon* 
   *force-float* *force-double*

   ;; fortran-atoms -- nothing is exported
   
   ;; numeric-vector
   
   numeric-vector-like numeric-vector elements nv-array-type 
   make-nv* make-nv create-nv copy-elements copy-nv-elements
   copy-nv

   ;; numeric-vector-wrappers -- nothing is exported

   ;; matrix-classes

   dense-matrix-like nrow ncol cm-index2 square-matrix-p square-matrix
   matrix-class matrix-kind restricted-elements set-restricted

   ;; compact-matrix

   dense-matrix upper-triangular-matrix lower-triangular-matrix hermitian-matrix

   ;; factorizations
   
   matrix-factorization component reconstruct lu-factorization lu-matrix ipiv
   qr-factorization qr-matrix cholesky-factorization factor hermitian-factorization

   ;; diagonal
   
   diagonal diagonal-integer diagonal-single diagonal-double
   diagonal-complex-single diagonal-complex-double make-diagonal
   create-diagonal nv->diagonal diagonal->nv matrix->diagonal
   diagonal->matrix

   ;; printing
   
   *pring-matrix-aligned* *print-matrix-padding*
   *print-matrix-precision*
   
   ;; matrix-operations
   
   submatrix make-matrix* make-matrix create-matrix reshape
   vector->column vector->row transpose stack-horizontally
   stack-vertically eye

   ;; readmacros

   read-vector-or-matrix v-syntax

   ;; fortran-call -- nothing is exported

   ;; linear-algebra

   mmm mm lu solve invert eigen least-squares least-squares-xx-inverse
   constrained-least-squares cholesky svd hermitian-factorization tr rank

   ;; adjustable

   *default-expansion* default-expansion adjustable size capacity add shrink
   adjustable-numeric-vector make-anv row-adjustable-matrix make-ra-matrix
   
   ))

;;;; ??? !!! maybe an lla-user package for playing around, once things
;;;; stabilize
