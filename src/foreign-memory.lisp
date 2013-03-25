;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla)

;;;; Lambda forms for accessing foreign memory.
;;;
;;; These are used to define the functions that access values in foreign memory.  It would be really convenient to use CFFI:MEM-AREF directly, but unfortunately C99 complex types are not (yet) supported (as of Oct 2011).  When they are available directly, we can eliminate these functions completely.  In the meantime, the compiler should be able to reduce the lambda forms, creating efficient code.

(eval-when (:compile-toplevel :load-toplevel :execute)
  (deftype maximum-array-size ()
    `(integer 0 #.(floor most-positive-fixnum (foreign-type-size :double))))
  (defun value-to-memory% (internal-type)
    "Return a (LAMBDA (POINTER INDEX VALUE) ...) form that can be used to write an element to an array in memory."
    `(lambda (pointer index value)
       (declare (type maximum-array-size index))
       ,(eswitch (internal-type)
          (+single+
           `(setf (mem-aref pointer :float index)
                  (coerce value 'single-float)))
          (+double+
           `(setf (mem-aref pointer :double index)
                  (coerce value 'double-float)))
          (+complex-single+
           `(let ((index2 (the fixnum (* 2 index))))
              (setf (mem-aref pointer :float index2)
                    (coerce (realpart value) 'single-float)
                    (mem-aref pointer :float (the fixnum (1+ index2)))
                    (coerce (imagpart value) 'single-float))))
          (+complex-double+
           `(let ((index2 (the fixnum (* 2 index))))
              (setf (mem-aref pointer :double index2)
                    (coerce (realpart value) 'double-float)
                    (mem-aref pointer :double (the fixnum (1+ index2)))
                    (coerce (imagpart value) 'double-float))))
          (+integer+
           `(setf (mem-aref pointer #-lla::int64 :int32
                                    #+lla::int64 :int64
                                    index)
                  (coerce value
                          '(signed-byte #-lla::int64 32 #+lla::int64 64)))))
       (values)))

  (defun value-from-memory% (internal-type)
    "Return a (LAMBDA (POINTER INDEX) ...) form that can be used to read an element from an array in memory."
    `(lambda (pointer index)
       (declare (type maximum-array-size index))
       ,(eswitch (internal-type)
          (+single+
           `(mem-aref pointer :float index))
          (+double+
           `(mem-aref pointer :double index))
          (+complex-single+
           `(let ((index2 (the fixnum (* 2 index))))
              (complex (mem-aref pointer :float index2)
                       (mem-aref pointer :float (1+ index2)))))
          (+complex-double+
           `(let ((index2 (the fixnum (* 2 index))))
              (complex (mem-aref pointer :double index2)
                       (mem-aref pointer :double (1+ index2)))))
          (+integer+
           `(mem-aref pointer #-lla::int64 :int32
                              #+lla::int64 :int64
                              index))))))

;;;; Accessing atoms individually.

(defun foreign-size (type)
  "Return the size of an internal type in bytes."
  (eswitch (type)
    (+integer+ (load-time-value (foreign-type-size #-lla::int64 :int32
                                                   #+lla::int64 :int64)))
    (+single+ (load-time-value (foreign-type-size :float)))
    (+double+ (load-time-value (foreign-type-size :double)))
    (+complex-single+ (load-time-value (* 2 (foreign-type-size :float))))
    (+complex-double+ (load-time-value (* 2 (foreign-type-size :double))))))

(defmacro define-foreign-aref ()
  `(progn
     (defun foreign-aref (pointer internal-type &optional (index 0))
       "Accessing 1d arrays in memory."
       (ecase internal-type
         ,@(loop for type in +internal-types+ collect
                 `(,type (,(value-from-memory% type) pointer index)))))
     (defun (setf foreign-aref) (value pointer internal-type
                                 &optional (index 0))
       "Setf expander for FOREIGN-AREF."
       (ecase internal-type
         ,@(loop for type in +internal-types+ collect
                 `(,type (,(value-to-memory% type) pointer index
                          (coerce value ',(lisp-type type)))))))))

(define-foreign-aref)

(defmacro with-fortran-atom ((pointer value type output) &body body)
  "Allocate memory for internal TYPE and set it to VALUE for body, which can use POINTER to access it.  When OUTPUT is given, the value is assigned to it after BODY.  The atom is automatically coerced to the correct type (by FOREIGN-AREF)."
  (check-type pointer symbol)
  (once-only (type)
    `(with-foreign-pointer (,pointer (foreign-size ,type))
       (setf (foreign-aref ,pointer ,type) ,value)
       ,@body
       ,@(when output
           `((setf ,output (foreign-aref ,pointer ,type)))))))

(defmacro with-fortran-character ((pointer character) &body body)
  "Make character available in an allocated memory area at POINTER for BODY."
  (check-type pointer symbol)
  (once-only (character)
    `(with-foreign-pointer (,pointer (foreign-type-size :char))
       (check-type ,character standard-char)
       (setf (mem-aref ,pointer :char) (char-code ,character))
       ,@body)))

;;;; Array copying functions
;;;
;;; Macros below define optimized copying for certain pairs of LLA types and Common Lisp array element types.  The functions below return the default ones for LLA.

;;;; Helper macros and standard specifications lists
;;;
;;; Note that when LLA copies from memory, it always aims for the same element type, so the other pairs need not be defined.

(eval-when (:compile-toplevel :load-toplevel :execute)

  (defun expand-specifications% (clause specifications)
    "Expand specifications using (clause internal-type element-type)."
    `(cond
       ,@(mapcan (lambda+ ((internal-type &rest element-types))
                   (mapcar (lambda (element-type)
                             (funcall clause internal-type element-type))
                           element-types))
                 specifications)
       (t (error 'lla-internal-error :message "Unhandled case."))))

  (defun all-to-specifications% ()
    "Return an optimization specification for all functions that copy to foreign memory."
    '((#.+single+ lla-single *)
      (#.+double+ lla-single lla-double *)
      (#.+complex-single+ lla-single lla-double lla-complex-single *)
      (#.+complex-double+ lla-single lla-double lla-complex-single
       lla-complex-double *)
      (#.+integer+ lla-integer *)))

  (defun all-from-specifications% ()
    "Return an optimization specification for all functions that copy from foreign memory."
    '((#.+single+ lla-single *)
      (#.+double+ lla-double *)
      (#.+complex-single+ lla-complex-single *)
      (#.+complex-double+ lla-complex-double *)
      (#.+integer+ lla-integer *)))

  (defmacro array-clause% ((array internal-type clause-element-type
                            clause-internal-type)
                           &body body)
    "Macro that generates a lambda form that can bed used in EXPAND-SPECIFICATIONS%."
    (with-gensyms (generic? array-type)
      `(lambda (,clause-internal-type ,clause-element-type)
         (let* ((,generic? (eq ,clause-element-type '*))
                (,array-type `(,(if ,generic? 'array 'simple-array)
                               ,,clause-element-type *)))
           `((and (typep ,',array ',,array-type)
                  (= ,',internal-type ,,clause-internal-type))
             (locally (declare ,@(unless ,generic?
                                   '((optimize speed (safety 0))))
                               (type ,,array-type ,',array))
               ,,@body)))))))

;;;; Copying to and from memory

(defun copy-array-to-memory (array pointer internal-type)
  "Copy the contents of ARRAY to the memory area of type INTERNAL-TYPE at POINTER."
  (declare (type internal-type internal-type)
           #+sbcl (sb-ext:muffle-conditions sb-ext:compiler-note))
  (check-type array array)
  (let ((size (array-total-size array)))
    (expanding
      (expand-specifications%
       (array-clause% (array internal-type element-type%
                             internal-type%)
         `(loop for index below size do
                   (,(value-to-memory% internal-type%) pointer index
                    (row-major-aref array index))))
       (all-to-specifications%))))
  (values))

(defun copy-array-from-memory (array pointer internal-type)
  "Copy the memory area of type INTERNAL-TYPE at POINTER to ARRAY."
  (declare (type internal-type internal-type))
  (check-type array array)
  (let+ ((size (array-total-size array)))
    (expanding
      (expand-specifications%
       (array-clause% (array internal-type element-type% internal-type%)
         `(loop for index below size do
           (setf (row-major-aref array index)
                 (coerce
                  (,(value-from-memory% internal-type%) pointer index)
                  ',(if (eq element-type% '*)
                        t
                        element-type%)))))
       (all-from-specifications%))))
  (values))

(defun create-array-from-memory (pointer internal-type dimensions
                                 &optional (element-type (lisp-type
                                                          internal-type)))
  "Create an array from contents at POINTER."
  (aprog1 (make-array dimensions :element-type element-type)
    (copy-array-from-memory it pointer internal-type)))

;;;; Transposing matrices to and from memory.

(defun transpose-matrix-to-memory (matrix pointer internal-type)
  "Transpose the contents of ARRAY to the memory area of type INTERNAL-TYPE at POINTER.  VECTORs are also handled."
  (etypecase matrix
    (vector (copy-array-to-memory matrix pointer internal-type))
    (aops:array-matrix
     (let+ (((nrow ncol) (array-dimensions matrix))
            (index 0))
       (declare (type fixnum index)
                #+sbcl (sb-ext:muffle-conditions sb-ext:compiler-note))
       (expanding
         (expand-specifications%
          (array-clause% (matrix internal-type element-type%
                                 internal-type%)
            `(loop for col-index fixnum below ncol do
                      (loop for row-index fixnum below nrow do
                               (,(value-to-memory% internal-type%) pointer index
                                (aref matrix row-index col-index))
                               (incf index))))
          (all-to-specifications%))))))
  (values))

(defun transpose-matrix-from-memory (matrix pointer internal-type)
  "Transpose the contents of ARRAY to the memory area of type INTERNAL-TYPE at POINTER.  VECTORs are also handled."
  (etypecase matrix
    (vector (copy-array-from-memory matrix pointer internal-type))
    (aops:array-matrix
     (let+ (((nrow ncol) (array-dimensions matrix))
            (index 0))
       (declare (type fixnum index))
       (expanding
         (expand-specifications%
          (array-clause% (matrix internal-type element-type%
                                 internal-type%)
            `(loop for col-index fixnum below ncol do
                      (loop for row-index fixnum below nrow do
                               (setf (aref matrix row-index col-index)
                                     (coerce
                                      (,(value-from-memory% internal-type%)
                                       pointer index)
                                      ',(if (eq element-type% '*)
                                            t
                                            element-type%)))
                               (incf index))))
          (all-from-specifications%))))))
  (values))

(defun create-transposed-matrix-from-memory (pointer internal-type dimensions
                                             &optional (element-type
                                                        (lisp-type
                                                         internal-type)))
  "Create a matrix from transposed contents at POINTER."
  (aprog1 (make-array dimensions :element-type element-type)
    (transpose-matrix-from-memory it pointer internal-type)))
