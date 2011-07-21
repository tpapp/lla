;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla)

(defparameter *lla-double?* t
  "Determines whether rational->float conversions result in double or single
   floats.")

(defun common-float-type (&rest objects)
  "Determine common float type for OBJECTS.  For use in LAPACK/BLAS calls."
  (common-lla-type objects :force-float? t :double? *lla-double?*))

(defun lookup-function (name)
  "Lookup pointer for NAME, asserting that it was found."
  (aprog1 (foreign-symbol-pointer name)
    (assert it () "~A not found." name)))

(defmacro lookup-procedure (lla-type name &optional (complex-name name))
  "Return an expression that returns a string for the LAPACK/BLAS procedure
name.  LLA-TYPE has to evaluate to a symbol denoting a float LLA type.  If you
need conditionals etc for the function name, do that outside this macro.  For
some functions (usually those involving Hermitian matrices), the names
actually differ based on whether the matrix is real or complex, use
COMPLEX-NAME in that case."
  ;; (check-type name symbol*)
  ;; (check-type complex-name symbol*)
  (flet ((lookup (letter name)
           `(load-time-value
             (lookup-function ,(format nil "~(~A~A_~)" letter name)))))
    `(ecase ,lla-type
       (:single ,(lookup "s" name))
       (:double ,(lookup "D" name))
       (:complex-single ,(lookup "C" complex-name))
       (:complex-double ,(lookup "Z" complex-name)))))

;;; interface

(defgeneric process-form (form environment)
  (:documentation "Return a list of argument specifications (atoms are
  converted into lists).")
  (:method (form environment)
    (macroexpand form environment )))

(defun process-forms (forms environment)
  "Process forms and return a list of argument specifications.  A form may
correspond to multiple arguments."
  (reduce #'append forms
          :key (lambda (f) (ensure-list (process-form f environment)))))

(defgeneric wrap-argument (argument pass parameters body)
  (:documentation "Return BODY wrapped in an environment generated for
  ARGUMENT in a given PASS.")
  (:method (argument pass parameters body)
    ;; default: just pass through body
    body))

(defun wrap-arguments (arguments pass parameters body)
  "Wrap BODY in arguments."
  (if arguments
      (wrap-argument (car arguments) pass parameters
                     (wrap-arguments (cdr arguments) pass parameters body))
      body))

(defgeneric argument-pointer (argument)
  (:documentation "Return the pointer for argument."))

(defun argument-pointers (arguments)
  "Return the list of pointers for all the arguments."
  (mapcar #'argument-pointer arguments))

(defun maybe-default-type (type parameters)
  "Return default type from parameters when TYPE is NIL."
  (aif type
       it
       (getf parameters :default-type)))

;;; implementation of specific types

;;; superclass of all pointers

(defstruct fortran-argument
  (pointer (gensym)))

(defmethod argument-pointer ((argument fortran-argument))
  (fortran-argument-pointer argument))

(defstruct (fortran-output (:include fortran-argument))
  (output nil))

(defmethod wrap-argument ((argument fortran-output) (pass (eql 'bindings))
                          parameters body)
  (let ((variable (fortran-output-output argument)))
    (if (and variable (not (keywordp variable)))
        `(let (,variable)
           ,body)
        body)))

;;; null pointer

(defmethod process-form ((form null) environment)
  (make-fortran-argument :pointer '(null-pointer)))

;;; atoms

(defstruct (fortran-atom (:include fortran-output))
  value (type nil) (coerce? nil))

(defmethod wrap-argument ((argument fortran-atom) (pass (eql 'main))
                          parameters body)
  (let+ (((&structure fortran-atom- pointer value type output coerce?)
          argument))
    `(with-fortran-atom (,pointer ,value
                         ,(maybe-default-type type parameters)
                         ,output ,coerce?)
       ,body)))

(defmacro &atom (value &key type output coerce?)
  (make-fortran-atom :value value :type type :output output :coerce? coerce?))

(defmacro &atom* (value &key type output &environment env)
  "Version of atom that coerces to the desired type."
  (process-form `(&atom ,value :type ,type :output ,output :coerce? t) env))

(defmacro &char (value &environment env)
  "Shorthand for character atoms."
  (process-form `(&atom ,value :type :char) env))

(defmacro &integer (value &key output &environment env)
  "Shorthand for integer atom."
  (process-form `(&atom ,value :type :integer :output ,output) env))

(defmacro &integers (&rest values &environment env)
  "Shorthand for integer atoms which are not modified."
  (loop for value in values
        collect (process-form `(&integer ,value) env)))

(defmethod process-form ((form (eql 0)) env)
  (process-form '(&atom* 0) env))

(defmethod process-form ((form (eql 1)) env)
  (process-form '(&atom* 1) env))

(defmethod process-form ((form character) env)
  (check-type form standard-char)
  (process-form `(&atom ,form :type :char) env))

;;; arrays

(defstruct (fortran-array (:include fortran-output))
  value (type nil) (transpose? nil) (output-dimensions nil)
  (output-transpose? nil))

(defun parse-array-output-specification (specification)
  "Parse an output specification (a list a single variable) and return it as
values."
  (if (eq specification :copy)
      :copy
      (let+ (((output &key dimensions transpose?)
              (aif specification
                   (ensure-list it)
                   (list nil))))
        (values output dimensions transpose?))))

(defmacro &array (value &key type transpose? output)
  "Fortran ARRAY.  When given, OUTPUT can be a specification of the
form (OUTPUT &KEY DIMENSIONS TRANSPOSE?)."
  (let+ (((&values o od ot) (parse-array-output-specification output)))
    (make-fortran-array :value value :type type :transpose? transpose?
                        :output o :output-dimensions od
                        :output-transpose? ot)))

(defmethod wrap-argument ((argument fortran-array) (pass (eql 'main))
                          parameters body)
  (let+ (((&structure-r/o fortran-array- pointer value type transpose? output
                          output-dimensions output-transpose?) argument))
    `(with-pinned-array (,pointer
                         ,value
                         ,(aif type
                               it
                               (getf parameters :default-type))
                         ,transpose?
                         ,output
                         ,output-dimensions
                         ,output-transpose?)
       ,body)))

;;; output arrays

(defstruct (fortran-output-array (:include fortran-output))
  dimensions (type nil) (transpose? nil))

(defmacro &output (output dimensions &key type transpose?)
  (make-fortran-output-array :output output :dimensions dimensions
                             :type type :transpose? transpose?))

(defmethod wrap-argument ((argument fortran-output-array) (pass (eql 'main))
                          parameters body)
  (let+ (((&structure-r/o fortran-output-array- pointer output dimensions type
                          transpose?) argument))
    `(with-array-output (,pointer
                         ,output
                         ,(maybe-default-type type parameters)
                         ,dimensions
                         ,transpose?)
       ,body)))

;;; work arrays

(defstruct (fortran-work-area (:include fortran-argument))
  (type nil) size)

(defmacro &work (size &optional type)
  (make-fortran-work-area :type type :size size))

(defmethod wrap-argument ((argument fortran-work-area) (pass (eql 'main))
                          parameters body)
  (let+ (((&structure-r/o fortran-work-area- pointer type size) argument))
    `(with-work-area (,pointer ,(maybe-default-type type parameters) ,size)
       ,body)))

;;; error handling

(define-condition lapack-error (error)
  ;; !! write method for formatting the error message
  ((lapack-procedure :initarg :lapack-procedure :type symbol*
		     :documentation "The name of the procedure."))
  (:documentation "The LAPACK procedure returned a nonzero info
  code."))

(define-condition lapack-invalid-argument (lapack-error)
  ((position :initarg :position :type fixnum
             :documentation "Position of the illegal argument"))
  (:documentation "An argument to a LAPACK procedure had an illegal
  value.  Generally, this indicates a bug in LLA and should not
  happen."))

(define-condition lapack-failure (lapack-error)
  ((info :initarg :info :type fixnum
          :documentation "INFO corresponding to error message."))
  (:documentation "Superclass of all LAPACK errors with a positive INFO"))

(define-condition lapack-singular-matrix (lapack-failure) ())

(define-condition lla-incompatible-dimensions (error) ())

(defstruct (lapack-info (:include fortran-argument))
  (variable (gensym))
  condition)

(defmacro &info (&optional (condition ''lapack-failure))
  (make-lapack-info :condition condition))

(define-symbol-macro &info (&info))

(defun lapack-info-wrap-argument (argument body)
  (let+ (((&structure-r/o lapack-info- pointer variable condition) argument))
    `(let (,variable)
       (with-fortran-atom (,pointer 0 :integer ,variable nil)
         ,body)
       (cond
         ((minusp ,variable) (error 'lapack-invalid-argument
                                    :position (- ,variable)))
         ((plusp ,variable) (error ',condition :info ,variable))))))

(defmethod wrap-argument ((argument lapack-info) (pass (eql 'call))
                          parameters body)
  (lapack-info-wrap-argument argument body))

(defmethod wrap-argument ((argument lapack-info) (pass (eql 'query))
                          parameters body)
  (lapack-info-wrap-argument argument body))

;;; work area query

(defstruct (lapack-work-query-area (:include fortran-argument))
  size
  type)

(defstruct (lapack-work-query-size (:include fortran-argument))
  size)

(defmacro &work-query (&optional type)
  (let ((size (gensym)))
    (list (make-lapack-work-query-area :size size :type type)
          (make-lapack-work-query-size :size size))))

(defmethod wrap-argument ((argument lapack-work-query-area)
                          (pass (eql 'bindings)) parameters body)
  (assert (getf parameters :query?) () "Call macro does not support queries.")
  `(let (,(lapack-work-query-area-size argument))
     ,body))

(defmethod wrap-argument ((argument lapack-work-query-area)
                          (pass (eql 'query)) parameters body)
  (let+ (((&structure-r/o lapack-work-query-area- pointer size type)
          argument))
    `(progn
       (with-fortran-atom (,pointer 0 ,(maybe-default-type type parameters)
                                    ,size t)
         ,body)
       (setf ,size (as-integer ,size)))))

(defmethod wrap-argument ((argument lapack-work-query-size)
                          (pass (eql 'query)) parameters body)
  (let+ (((&structure-r/o lapack-work-query-size- pointer) argument))
    `(with-fortran-atom (,pointer -1 :integer nil nil)
       ,body)))

(defmethod wrap-argument ((argument lapack-work-query-area)
                          (pass (eql 'call)) parameters body)
  (let+ (((&structure-r/o lapack-work-query-area- pointer size type)
          argument))
    `(with-work-area (,pointer ,(maybe-default-type type parameters) ,size)
       ,body)))

(defmethod wrap-argument ((argument lapack-work-query-size)
                          (pass (eql 'call)) parameters body)
  (let+ (((&structure-r/o lapack-work-query-size- pointer size) argument))
    `(with-fortran-atom (,pointer ,size :integer nil nil)
       ,body)))

;;; various call interfaces

(defun function-names (name)
  "Parse LAPACK/BLAC fuction names, return a list of two strings (real and
complex)."
  (let+ (((name &optional (complex-name name)) (ensure-list name)))
    (check-types (name complex-name) string)
    (list name complex-name)))

(defun fortran-call-form (function-pointer arguments)
  "Return a form that calls FUNCTION-POINTER with given arguments."
  `(foreign-funcall-pointer ,function-pointer ()
                            ,@(loop for arg in arguments
                                    appending
                                    `(:pointer ,(argument-pointer arg)))
                            :void))

(defun named-call-form (name type-var arguments)
  (fortran-call-form `(lookup-procedure ,type-var ,@(function-names name))
                     arguments))

(defmacro blas-call ((name type value) &body forms &environment env)
  "!!"
  (let* ((type-var (gensym "TYPE"))
         (arguments (process-forms forms env))
         (parameters `(:default-type ,type-var)))
    `(let ((,type-var ,type))
       ,(wrap-arguments arguments 'bindings parameters
                        `(progn
                           ,(wrap-arguments arguments 'main parameters
                                            (named-call-form name type-var
                                                             arguments))
                           ,value)))))

(defmacro lapack-call ((name type value) &body forms &environment env)
  "!!"
  (let* ((type-var (gensym "TYPE"))
         (arguments (process-forms forms env))
         (parameters `(:default-type ,type-var)))
    (assert (<= (count-if #'lapack-info-p arguments) 1))
    `(let ((,type-var ,type))
       ,(wrap-arguments
         arguments 'bindings parameters
         `(progn
            ,(wrap-arguments
              arguments 'main parameters
              (wrap-arguments
               arguments 'call parameters
               (named-call-form name type-var
                                arguments)))
            ,value)))))

(defmacro lapack-call-w/query ((name type value) &body forms &environment env)
  "!!"
  (let* ((type-var (gensym "TYPE"))
         (arguments (process-forms forms env))
         (parameters `(:default-type ,type-var :query? t))
         (call-form (named-call-form name type-var arguments)))
    (assert (<= (count-if #'lapack-info-p arguments) 1))
    `(let ((,type-var ,type))
       ,(wrap-arguments
         arguments 'bindings parameters
         `(progn
            ,(wrap-arguments
              arguments 'main parameters
              `(progn
                 ,(wrap-arguments arguments 'query parameters call-form)
                 ,(wrap-arguments arguments 'call parameters call-form)))
            ,value)))))

;;; floating point traps
;;;
;;; Apparently, the only trap that we need to mask is division by zero, and
;;; that only for a few operations.  Non-numerical floating points values are
;;; used internally (eg in SVD calculations), but only reals are returned.

#-(or sbcl cmu)
(defmacro with-lapack-traps-masked (&body body)
  (warn "No with-lapack-traps-masked macro provided for your ~
  implementation -- some operations may signal an error.")
  `(progn
     ,@body))

#+sbcl
(defmacro with-lapack-traps-masked (&body body)
  `(sb-int:with-float-traps-masked (:divide-by-zero :invalid)
     ,@body))

#+cmu
(defmacro with-lapack-traps-masked (&body body)
  `(extensions:with-float-traps-masked (:divide-by-zero :invalid)
     ,@body))
