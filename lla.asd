(defpackage #:lla-asd
  (:use :cl :asdf))

(in-package #:lla-asd)

(defparameter *fasl-directory*
  (make-pathname :directory '(:relative
			      #+sbcl "fasl-sbcl"
			      #+openmcl "fasl-ccl"
			      #+cmu "fasl-cmucl"
			      #+clisp "fasl-clisp"
			      #-(or sbcl openmcl clisp cmucl) "fasl"
			      )))

(defsystem #:lla
  :description "Lisp Linear Algebra"
  :version "alpha"
  :author "Tamas K Papp <tkpapp@gmail.com>"
  :license "BSD without advertising clause"
  :serial t
  :components 
  ((:module 
    "package-init"
    :pathname #P "src/"
    :components
    ((:file "package")))
   (:module 
    "basics"
    :pathname #P"src/"
    :depends-on ("package-init")
    :serial t
    :components
    ((:file "load-libs")
     (:file "blas-cffi")
     (:file "lapack-cffi")
     (:file "utilities")
     (:file "types")
     (:file "scalar")
     (:file "numeric-vector")
     (:file "matrix")
     (:file "fortran-call")
;; !! include these when they stabilize
     )))
  :depends-on
  (:cl-utilities :iterate :metabang-bind :cffi :cl-match :xarray
		 :lift))


;;;; !! ASDF loading for unit tests. 
;;;;
;;;; ?? I think it should go into a separate defsystem so it could be
;;;;    loaded separately -- Tamas
;;;;
;;;; ?? providing something else besides ASDF?  Mudballs? -- Tamas
