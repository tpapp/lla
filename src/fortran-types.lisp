(in-package :lla)

;; type conversions between Rif's fortran types and pointers

(eval-when (:compile-toplevel :load-toplevel :execute)
  (defctype fortran-logical :pointer)
  (defctype fortran-int :pointer)
  (defctype fortran-float :pointer)
  (defctype fortran-double :pointer)
  (defctype fortran-complex-float :pointer)
  (defctype fortran-complex-double :pointer)
  (defctype cffi-fnv-float :pointer)
  (defctype cffi-fnv-double :pointer)
  (defctype cffi-fnv-complex-float :pointer)
  (defctype cffi-fnv-complex-double :pointer)
  (defctype cffi-fnv-int32 :pointer))
