(in-package :lla)

;;;;
;;;;  This might require some tweaking, not (yet) tested on OSs other
;;;;  than Linux.
;;;;
;;;;  !! do I need gfortran?
;;;;  !! do I need an eval-when wrapper?  CFFI seems to take care of things neatly

(define-foreign-library :libblas
  (cffi-features:darwin   "libblas.dylib")
  (cffi-features:unix     (:or "libblas.so" "libblas.so.3gf"))
  (cffi-features:windows  "libblas.dll"))

(define-foreign-library :liblapack
  (cffi-features:darwin   "liblapack.dylib")
  (cffi-features:unix     (:or "liblapack.so" "liblapack.so.3gf" "liblapack_atlas.so.3gf"))
  (cffi-features:windows  "liblapack.dll"))

(load-foreign-library :libblas)
(load-foreign-library :liblapack)
