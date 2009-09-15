(in-package :lla)

;;;; ** direct memory access
;;;;
;;;; mem-aref* is similar to cffi:mem-aref, except for two things: (1)
;;;; type is lla-type, and this (2) it handles complex types, too.

;;;; ?? Currently, I did not really concentrate on optimizing these
;;;; things, as they will be used for occasional scalar access only.

(defun mem-aref* (ptr lla-type &optional (index 0))
  (check-type lla-type lla-type)
;;  (check-type index fixnum)
  "Accessing 1d arrays in memory.  Type is LLA-type."
  (ecase lla-type
    (:integer (mem-aref ptr :uint32 index))
    (:single (mem-aref ptr :float index))
    (:double (mem-aref ptr :double index))
    (:complex-single
       (complex (mem-aref ptr :float (* 2 index))
		(mem-aref ptr :float (1+ (* 2 index)))))
    (:complex-double
       (complex (mem-aref ptr :double (* 2 index))
		(mem-aref ptr :double (1+ (* 2 index)))))))

(defun (setf mem-aref*) (value ptr lla-type &optional (index 0))
  (check-type lla-type lla-type)
;;  (check-type index fixnum)
  "Setf expander for accessing 1d arrays in memory.  Type is LLA-type."
  (ecase lla-type
    (:integer (setf (mem-aref ptr :uint32 index) value))
    (:single (setf (mem-aref ptr :float index) value))
    (:double (setf (mem-aref ptr :double index) value))
    (:complex-single
       (setf (mem-aref ptr :float (* 2 index)) (realpart value)
	     (mem-aref ptr :float (1+ (* 2 index))) (imagpart value)))
    (:complex-double
       (setf (mem-aref ptr :double (* 2 index)) (realpart value)
	     (mem-aref ptr :double (1+ (* 2 index))) (imagpart value)))))

(defun foreign-size* (lla-type)
  "Return the size of a foreign type, in bytes."
  (check-type lla-type lla-type)
  (ecase lla-type
    (:integer #.(foreign-type-size :uint32))
    (:single #.(foreign-type-size :float))
    (:double #.(foreign-type-size :double))
    (:complex-single #.(* 2 (foreign-type-size :float)))
    (:complex-double #.(* 2 (foreign-type-size :double)))))
    

;;;; ** interfacing with Fortran
;;;;
;;;;

(defmacro with-fortran-scalar ((value pointer lla-type) &body body)
  "Allocate memory for lla-type and set it to value for body."
  (check-type pointer symbol)
  (once-only (value lla-type)
    `(with-foreign-pointer (,pointer (foreign-size* ,lla-type))
       (setf (mem-aref* ,pointer ,lla-type) ,value)
       ,@body)))

(with-multiple-bindings with-fortran-scalar)

;; ?? naming convention: above could be with-fortran-scalar-input, if
;; there were -input-output and -input version, but there won't be:
;; work area queries will be another specialized macro, and returning
;; a single calculated value (eg a matrix norm) will be another.
;; !! need to write those -- Tamas


