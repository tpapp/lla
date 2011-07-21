;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla)

(defun foreign-size* (type)
  "Return the size of an LLA type or :CHAR, in bytes."
  (ecase type
    (:char (load-time-value (foreign-type-size :char)))
    (:integer (load-time-value (foreign-type-size #-lla::int64 :int32
                                                  #+lla::int64 :int64)))
    (:single (load-time-value (foreign-type-size :float)))
    (:double (load-time-value (foreign-type-size :double)))
    (:complex-single (load-time-value (* 2 (foreign-type-size :float))))
    (:complex-double (load-time-value (* 2 (foreign-type-size :double))))))

(defun mem-aref* (ptr type &optional (index 0))
  "Accessing 1d arrays in memory.  TYPE is an LLA type or :CHAR."
  (ecase type
    (:char (code-char (mem-aref ptr :char index)))
    (:integer (mem-aref ptr #-lla::int64 :int32 #+lla::int64 :int64 index))
    (:single (mem-aref ptr :float index))
    (:double (mem-aref ptr :double index))
    (:complex-single
       (complex (mem-aref ptr :float (* 2 index))
		(mem-aref ptr :float (1+ (* 2 index)))))
    (:complex-double
       (complex (mem-aref ptr :double (* 2 index))
		(mem-aref ptr :double (1+ (* 2 index)))))))

(defun (setf mem-aref*) (value ptr type &optional (index 0))
  "Setf expander for accessing 1d arrays in memory.  TYPE is an LLA type
or :CHAR."
  (ecase type
    (:char (setf (mem-aref ptr :char index) (char-code value)))
    (:integer (setf (mem-aref ptr #-lla::int64 :int32
                              #+lla::int64 :int64 index) value))
    (:single (setf (mem-aref ptr :float index) value))
    (:double (setf (mem-aref ptr :double index) value))
    (:complex-single
       (setf (mem-aref ptr :float (* 2 index)) (realpart value)
	     (mem-aref ptr :float (1+ (* 2 index))) (imagpart value)))
    (:complex-double
       (setf (mem-aref ptr :double (* 2 index)) (realpart value)
	     (mem-aref ptr :double (1+ (* 2 index))) (imagpart value)))))

(defmacro with-fortran-atom ((pointer value type output coerce?) &body body)
  "Allocate memory for TYPE (types handled by FOREIGN-SIZE* etc) and set it to
VALUE for body, which can use POINTER to access it.  When OUTPUT is given, the
value is assigned to it after BODY.  When coerce?, the atom is coerced to the
correct type."
  (check-type pointer symbol)
  (once-only (type)
    `(with-foreign-pointer (,pointer (foreign-size* ,type))
       (setf (mem-aref* ,pointer ,type) ,(if coerce?
                                             `(coerce* ,value ,type)
                                             value))
       ,@body
       ,@(when output
           `((setf ,output (mem-aref* ,pointer ,type)))))))
