;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla)

;;; These macros take care of array pinning and the conversion of constants
;;; to pointers to an allocated memory area with the value so that everything
;;; is set up to call Fortran/BLAS/LAPACK functions.
;;;
;;; Since LAPACK functions are sometimes called two times (eg to query work
;;; area sizes), expansions take care in 'passes', such as
;;; 
;;;  - bindings: establishes the bindings (also empty variables)
;;;
;;;  - main: for arguments that are the same regardless of what kind of call
;;;    is made
;;; 
;;;  - query: for querying work area sizes
;;; 
;;;  - call: the actual function call
;;; 
;;; The DSL is implemented via macros which expand to structures, which are in
;;; then handed to WRAP-ARGUMENT for each pass.

;;;; generic interface and helper functions

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
  "Wrap BODY in arguments.  Convenienve function used to implement the
expansion."
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

;;;; implementation of specific types

;;; fortran-argument

(defstruct fortran-argument
  "Superclass of all arguments with pointers."
  (pointer (gensym)))

(defmethod argument-pointer ((argument fortran-argument))
  (fortran-argument-pointer argument))

;;; fortran-output

(defstruct (fortran-output (:include fortran-argument))
  "For parameters that have an output."
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

;;; characters

(defstruct (fortran-character (:include fortran-argument))
  "Characters passed to FORTRAN.  Input only, for specifying triangle
orientation, etc."
  value)

(defmethod wrap-argument ((argument fortran-character) (pass (eql 'main))
                          parameters body)
  (let+ (((&structure-r/o fortran-character- pointer value) argument))
    `(with-fortran-character (,pointer ,value)
       ,body)))

(defmacro &char (value)
  "Shorthand for character atoms."
  (make-fortran-character :value value))

(defmethod process-form ((form character) env)
  (check-type form standard-char)
  (process-form `(&char ,form) env))

;;; atoms

(defstruct (fortran-atom (:include fortran-output))
  "Atoms passed to FORTRAN."
  value (type nil) (coerce? nil))

(defmethod wrap-argument ((argument fortran-atom) (pass (eql 'main))
                          parameters body)
  (let+ (((&structure-r/o fortran-atom- pointer value type output)
          argument))
    `(with-fortran-atom (,pointer ,value ,(maybe-default-type type parameters)
                                  ,output)
       ,body)))

(defmacro &atom (value &key type output)
  "Atoms passed to FORTRAN.  When not given, TYPE is inferred from the call's
default.  VALUE is coerced to the desired type.  When OUTPUT is given, value
is read after the call and placed there."
  (make-fortran-atom :value value :type type :output output))

(defmacro &integer (value &key output &environment env)
  "Shorthand for integer atom."
  (process-form `(&atom ,value :type +integer+ :output ,output) env))

(defmacro &integers (&rest values &environment env)
  "Shorthand for integer atoms which are not modified."
  (loop for value in values
        collect (process-form `(&integer ,value) env)))

(defmethod process-form ((form (eql 0)) env)
  (process-form '(&atom 0) env))

(defmethod process-form ((form (eql 1)) env)
  (process-form '(&atom 1) env))


;;; arrays

(defstruct (fortran-array (:include fortran-output))
  "Arrays which are pinned."
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
  "Output arrays."
  dimensions (type nil) (transpose? nil))

(defmacro &output (output dimensions &key type transpose?)
  "A memory area of DIMENSIONS is allocated, and the contents are assigned to
OUTPUT after the call.  When not given, TYPE is the call's default.  When
TRANSPOSE?, the result is transposed."
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
  "Work area."
  (type nil) size)

(defmacro &work (size &optional type)
  "Allocate a work area of SIZE.  When TYPE is not given, the call's default
is used."
  (make-fortran-work-area :type type :size size))

(defmethod wrap-argument ((argument fortran-work-area) (pass (eql 'main))
                          parameters body)
  (let+ (((&structure-r/o fortran-work-area- pointer type size) argument))
    `(with-work-area (,pointer ,(maybe-default-type type parameters) ,size)
       ,body)))


(defstruct (lapack-info (:include fortran-argument))
  (variable (gensym))
  condition)

;;; call info

(defmacro &info (&optional (condition ''lapack-failure))
  "Argument for checking whether the call was executed without an error.
Automatically takes care of raising the appropriate condition if it wasn't.
CONDITION specifies the condition to raise in case of positive error codes."
  (make-lapack-info :condition condition))

(define-symbol-macro &info (&info))

(defun lapack-info-wrap-argument (argument body)
  (let+ (((&structure-r/o lapack-info- pointer variable condition) argument))
    `(let (,variable)
       (with-fortran-atom (,pointer 0 +integer+ ,variable)
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
;;; 
;;; &work-query expands to TWO structures which share a SIZE argument, they
;;; cooperate for the query.

(defstruct (lapack-work-query-area (:include fortran-argument))
  size
  type)

(defstruct (lapack-work-query-size (:include fortran-argument))
  size)

(defmacro &work-query (&optional type)
  "Work area query, takes the place of TWO fortran arguments."
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
                                    ,size)
         ,body)
       (setf ,size (as-integer ,size)))))

(defmethod wrap-argument ((argument lapack-work-query-size)
                          (pass (eql 'query)) parameters body)
  (let+ (((&structure-r/o lapack-work-query-size- pointer) argument))
    `(with-fortran-atom (,pointer -1 +integer+ nil)
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
    `(with-fortran-atom (,pointer ,size +integer+ nil)
       ,body)))

;;; various call interfaces

(defun blas-lapack-function-name (type name)
  "Return the BLAS/LAPACK foreign function name.  TYPE is the internal type,
NAME is one of the following: NAME, (NAME), which are used for both complex
and real names, or (REAL-NAME COMPLEX-NAME)."
  (let+ (((real-name &optional (complex-name name)) (ensure-list name))
         (letter (switch (type)
                   (+single+ "S")
                   (+double+ "D")
                   (+complex-single+ "C")
                   (+complex-double+ "Z")))
         (name (if (complex? type)
                   complex-name
                   real-name)))
    (format nil "~(~A~A_~)" letter name)))

(defun arguments-for-cffi (arguments)
  "Return a list that can be use in a CFFI call."
  (loop for arg in arguments appending `(:pointer ,(argument-pointer arg))))

(defun blas-lapack-call-form (type-var name arguments)
  "Return a form BLAS/LAPACK calls, conditioning on TYPE-VAR.  See
BLAS-LAPACK-FUNCTION-NAME for the interpretation of "
  (let ((arguments (arguments-for-cffi arguments)))
    `(ecase ,type-var
       ,@(loop for type in +float-types+
               collect
               `(,type
                 (cffi:foreign-funcall 
                  ,(blas-lapack-function-name type name)
                  ,@arguments
                  :void))))))

;;;; Main interface
;;; 
;;; Common conventions:
;;;
;;;  1. NAME is either a string or a list of two strings (real/complex)
;;; 
;;;  2. VALUE is the form returned after the call

(defmacro blas-call ((name type value) &body forms &environment env)
  "BLAS call."
  (let* ((type-var (gensym "TYPE"))
         (arguments (process-forms forms env))
         (parameters `(:default-type ,type-var)))
    `(let ((,type-var ,type))
       ,(wrap-arguments
         arguments 'bindings parameters
         `(progn
            ,(wrap-arguments arguments 'main parameters
                             (blas-lapack-call-form type-var name arguments))
            ,value)))))

(defmacro lapack-call ((name type value) &body forms &environment env)
  "LAPACK call, takes an &info argument."
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
              (wrap-arguments arguments 'call parameters
                              (blas-lapack-call-form type-var name
                                                     arguments)))
            ,value)))))

(defmacro lapack-call-w/query ((name type value) &body forms &environment env)
  "LAPACK call which also takes &work-query arguments (in place of two FORTRAN
arguments)."
  (let* ((type-var (gensym "TYPE"))
         (arguments (process-forms forms env))
         (parameters `(:default-type ,type-var :query? t))
         (call-form (blas-lapack-call-form type-var name arguments)))
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

;;;; floating point traps
;;;
;;; Apparently, the only trap that we need to mask is division by zero, and
;;; that only for a few operations.  Non-numerical floating points values are
;;; used internally (eg in SVD calculations), but only reals are returned.

#-(or sbcl cmu)
(defmacro with-fp-traps-masked (&body body)
  (warn "No with-lapack-traps-masked macro provided for your ~
  implementation -- some operations may signal an error.")
  `(progn
     ,@body))

#+sbcl
(defmacro with-fp-traps-masked (&body body)
  `(sb-int:with-float-traps-masked (:divide-by-zero :invalid)
     ,@body))

#+cmu
(defmacro with-fp-traps-masked (&body body)
  `(extensions:with-float-traps-masked (:divide-by-zero :invalid)
     ,@body))
