(in-package :lla)

;;; Some character constants as integers.  In theory, could be defined using
;;; CHAR-CODE on ASCII systems.

(defconstant +c+ 67 "Numerical code for character C (ASCII), for use in LAPACK.")
(defconstant +l+ 76 "Numerical code for character L (ASCII), for use in LAPACK.")
(defconstant +n+ 78 "Numerical code for character N (ASCII), for use in LAPACK.")
(defconstant +t+ 84 "Numerical code for character T (ASCII), for use in LAPACK.")
(defconstant +u+ 85 "Numerical code for character U (ASCII), for use in LAPACK.")

;;;; ** direct memory access
;;;;
;;;; mem-aref* is similar to cffi:mem-aref, except for two things: (1)
;;;; type an LLA type or :char, and (2) it handles complex types, too.

;;;; ?? Currently, I did not really concentrate on optimizing these
;;;; things, as they will be used for occasional scalar access only.
;;;;
;;;; !! Once CFFI gets to handle complex types, life will be so much
;;;; !! simpler and we can dispense with a lot of this code.

(defun foreign-size* (lla-type)
  "Return the size of an LLA-TYPE or :CHAR, in bytes."
  (ecase lla-type
    (:char (load-time-value (foreign-type-size :char)))
    (:integer (load-time-value (foreign-type-size #-lla-int64 :uint32
                                                  #+lla-int64 :uint64)))
    (:single (load-time-value (foreign-type-size :float)))
    (:double (load-time-value (foreign-type-size :double)))
    (:complex-single (load-time-value (* 2 (foreign-type-size :float))))
    (:complex-double (load-time-value (* 2 (foreign-type-size :double))))))

(defun mem-aref* (ptr lla-type &optional (index 0))
  "Accessing 1d arrays in memory.  TYPE is LLA-TYPE or :CHAR."
  (ecase lla-type
    (:char (code-char (mem-aref ptr :char index)))
    (:integer (mem-aref ptr #-lla-int64 :uint32 #+lla-int64 :uint64 index))
    (:single (mem-aref ptr :float index))
    (:double (mem-aref ptr :double index))
    (:complex-single
       (complex (mem-aref ptr :float (* 2 index))
		(mem-aref ptr :float (1+ (* 2 index)))))
    (:complex-double
       (complex (mem-aref ptr :double (* 2 index))
		(mem-aref ptr :double (1+ (* 2 index)))))))

(defun (setf mem-aref*) (value ptr type &optional (index 0))
  "Setf expander for accessing 1d arrays in memory.  TYPE is LLA-TYPE
or :CHAR"
  (ecase type
    (:char (setf (mem-aref ptr :char index) (char-code value)))
    (:integer (setf (mem-aref ptr #-lla-int64 :uint32
                              #+lla-int64 :uint64 index) value))
    (:single (setf (mem-aref ptr :float index) value))
    (:double (setf (mem-aref ptr :double index) value))
    (:complex-single
       (setf (mem-aref ptr :float (* 2 index)) (realpart value)
	     (mem-aref ptr :float (1+ (* 2 index))) (imagpart value)))
    (:complex-double
       (setf (mem-aref ptr :double (* 2 index)) (realpart value)
	     (mem-aref ptr :double (1+ (* 2 index))) (imagpart value)))))

;; (defmacro with-foreign-atom ((pointer type value) &body body)
;;   "Allocate memory for TYPE (types handled by FOREIGN-SIZE* etc) and
;; set it to VALUE for body, which can use POINTER to access it. VALUE
;; and POINTER can be the same symbol, this is the default if VALUE is
;; omitted."
;;   (check-type pointer symbol)
;;   (once-only (type)
;;     (with-unique-names (value-saved)
;;       `(let ((,value-saved ,value)) ; value can be the same as pointer
;;          (with-foreign-pointer (,pointer (foreign-size* ,type))
;;            (setf (mem-aref* ,pointer ,type) ,value-saved)
;;            ,@body)))))

;; (define-with-multiple-bindings with-foreign-atom)
