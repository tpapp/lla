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
    :pathname #P"src/"
    :components
    ((:file "package")))
   (:module
    "fortran-interface"
    :pathname #P"src/"
    :serial t
    :components
    ((:file "architecture")
     (:file "load-libs")
     (:file "fortran-types")
     (:file "blas-cffi")
     (:file "lapack-cffi")))
   (:module 
    "basics"
    :pathname #P"src/"
    :depends-on ("package-init" "fortran-interface")
    :serial t
    :components
    ((:file "utilities")
     (:file "types")
     (:file "fortran-atoms")
     (:file "numeric-vector")
     (:file "numeric-vector-wrappers")
     (:file "matrix-classes")
     (:file "compact-matrix")
     (:file "matrix-operations")
     (:file "factorizations")
     (:file "submatrix")
     (:file "bind-extensions")
     (:file "diagonal")
     (:file "printing")
     (:file "readmacros")
     (:file "creation-conversion")
     (:file "binary-operations")))
   (:module
    "linear-algebra"
    :pathname #P"src/"
    :depends-on ("basics")
    :serial t
    :components
     ((:file "fortran-call")
      (:file "linear-algebra")))
   (:module
    "extensions"
    :pathname #P"src/"
    :depends-on ("basics")
    :serial t
    :components
    ((:file "adjustable"))))
  :depends-on
  (:cl-utilities :iterate :metabang-bind :cffi :xarray
                 :anaphora :named-readtables :alexandria))

;;;; ?? providing something else besides ASDF?  Mudballs? -- Tamas
