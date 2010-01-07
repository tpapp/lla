(in-package :lla)

;;;; ** direct memory access
;;;;
;;;; mem-aref* is similar to cffi:mem-aref, except for two things: (1)
;;;; type is lla-type or :char, and (2) it handles complex types, too.

;;;; ?? Currently, I did not really concentrate on optimizing these
;;;; things, as they will be used for occasional scalar access only.
;;;;
;;;; !! Once CFFI gets to handle complex types, life will be so much
;;;; !! simpler and we can dispense with a lot of this code.

(defun foreign-size* (type)
  "Return the size of an LLA-TYPE or :CHAR, in bytes."
  (ecase type
    (:char (load-time-value (foreign-type-size :char)))
    (:integer (load-time-value (foreign-type-size :uint32)))
    (:single (load-time-value (foreign-type-size :float)))
    (:double (load-time-value (foreign-type-size :double)))
    (:complex-single (load-time-value (* 2 (foreign-type-size :float))))
    (:complex-double (load-time-value (* 2 (foreign-type-size :double))))))

(defun mem-aref* (ptr type &optional (index 0))
;;  (check-type index fixnum)
  "Accessing 1d arrays in memory.  TYPE is LLA-TYPE or :CHAR."
  (ecase type
    (:char (code-char (mem-aref ptr :char index)))
    (:integer (mem-aref ptr :uint32 index))
    (:single (mem-aref ptr :float index))
    (:double (mem-aref ptr :double index))
    (:complex-single
       (complex (mem-aref ptr :float (* 2 index))
		(mem-aref ptr :float (1+ (* 2 index)))))
    (:complex-double
       (complex (mem-aref ptr :double (* 2 index))
		(mem-aref ptr :double (1+ (* 2 index)))))))

(defun (setf mem-aref*) (value ptr type &optional (index 0))
;;  (check-type index fixnum)
  "Setf expander for accessing 1d arrays in memory.  TYPE is LLA-TYPE or :CHAR"
  (ecase type
    (:char (setf (mem-aref ptr :char index) (char-code value)))
    (:integer (setf (mem-aref ptr :uint32 index) value))
    (:single (setf (mem-aref ptr :float index) value))
    (:double (setf (mem-aref ptr :double index) value))
    (:complex-single
       (setf (mem-aref ptr :float (* 2 index)) (realpart value)
	     (mem-aref ptr :float (1+ (* 2 index))) (imagpart value)))
    (:complex-double
       (setf (mem-aref ptr :double (* 2 index)) (realpart value)
	     (mem-aref ptr :double (1+ (* 2 index))) (imagpart value)))))


;;;;  interfacing with Fortran
;;;;
;;;;

(defmacro with-fortran-atom ((type pointer &optional (value pointer)) &body body)
  "Allocate memory for TYPE (types handled by FOREIGN-SIZE* etc) and
set it to VALUE for body, which can use POINTER to access it.  VALUE
and POINTER can be the same symbol, this is the default if VALUE is
omitted."
  (check-type pointer symbol)
  (once-only (type)
    (with-unique-names (value-saved)
      `(let ((,value-saved ,value))     ; value can be the same as pointer
         (with-foreign-pointer (,pointer (foreign-size* ,type))
           (setf (mem-aref* ,pointer ,type) ,value-saved)
           ,@body)))))

(with-multiple-bindings with-fortran-atom)

;; ?? naming convention: above could be with-fortran-scalar-input, if
;; there were -input-output and -input version, but there won't be:
;; work area queries will be another specialized macro, and returning
;; a single calculated value (eg a matrix norm) will be another.
;; !! need to write those -- Tamas


;;;; this is missing from CFFI at the moment
(with-multiple-bindings with-foreign-pointer)
