(in-package #:lla-asd)

(defpackage #:lla
  (:use :common-lisp
	:cl-utilities
        :iterate
        :bind
	:cffi
	:cl-match
	:xarray)
  (:shadowing-import-from :iterate :collecting :collect)
  (:export 
   ;; !! exported interface should/will be refined and changed,
   ;; !! currently it is just a sketch

   ;; types
   
   ;; scalar
   
   mem-aref* foreign-size* with-fortran-scalar with-fortran-scalars
   
   ;; numeric-vector
   
   numeric-vector numeric-vector-data numeric-vector-integer
   numeric-vector-single numeric-vector-double
   numeric-vector-complex-single numeric-vector-complex-double
   nv-element-type make-numeric-vector nv-copy element-lla-type

   ;; matrix
   
   matrix nrow ncol data *pring-matrix-aligned* *print-matrix-padding*
   *print-matrix-precision* dense-matrix make-matrix vector->matrix-col
   vector->matrix-row matrix->vector lu ipiv
   
   ;; fortran-call
   
   lu solve mm eigen

   ))

(defpackage #:lla-unit-tests
  (:documentation "Unit, validation, and regression testing for lla")
  (:use :common-lisp
	:lift
	:lla)
  (:export run-lla-tests))

;;;; ??? !!! maybe an lla-user package for playing around, once things
;;;; stabilize
