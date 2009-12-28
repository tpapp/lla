(in-package :lla)

;;;;  Code for interfacing with Fortran calls.


;;;;  Helper functions that generate the correct LAPACK/BLAS function
;;;;  names based on a "root" function name.  For some functions
;;;;  (usually those involving Hermitian matrices), the roots actually
;;;;  differ based on whether the matrix is real or complex.

(declaim (inline lb-procedure-name lb-procedure-name2))
(defun lb-procedure-name (name lla-type)
  "Evaluate to the LAPACK/BLAS procedure name.  LLA-TYPE has to
evaluate to a symbol.  If you need conditionals etc, do that outside
this macro."
  (check-type name symbol)
  (check-type lla-type symbol)
  (ecase lla-type
    (:single (make-symbol* "%S" name))
    (:double (make-symbol* "%D" name))
    (:complex-single (make-symbol* "%C" name))
    (:complex-double (make-symbol* "%Z" name))))

(defun lb-procedure-name2 (real-name complex-name lla-type)
  "Evaluate to the LAPACK/BLAS procedure name, differentiating real
and complex cases.  :single or :double is returned as the second value
as appropriate, and the third value is true iff lla-type is complex."
  (check-type real-name symbol)
  (check-type complex-name symbol)
  (check-type lla-type symbol)
  (ecase lla-type
    (:single (values (make-symbol* "%S" real-name) :single nil))
    (:double (values (make-symbol* "%D" real-name) :double nil))
    (:complex-single (values (make-symbol* "%C" complex-name) :single t))
    (:complex-double (values (make-symbol* "%Z" complex-name) :double t))))



;;;;  Some LAPACK procedures can signal errors, they do this via an
;;;;  INFO output integer.  Here we provide macros to capture these
;;;;  errors.
;;;;
;;;;  Currently, a general lapack-error condition is thrown.  This may
;;;;  be too coarse, as INFO actually provides quite a bit of
;;;;  information of what went wrong.  There are usually two kinds of
;;;;  errors: invalid/nonsensical arguments (should never happen), and
;;;;  incomputable problems (matrix being singular, usually with
;;;;  detailed info).  !!! We should capture the latter and provide
;;;;  more sensible error messages.  Maybe instead of throwing
;;;;  LAPACK-ERROR, macros should accept a form that tells what kind
;;;;  of error to signal.  Expect CALL-WITH-INFO-CHECK to be modified
;;;;  in the future, eg all current arguments go into a list, then
;;;;  body gets PROCEDURE and INFO values in case of an error and may
;;;;  do what it wants with them.
;;;;
;;;;  ??? Most (maybe all?) LAPACK functions just have INFO as their
;;;;  last argument, so WITH-INFO-CHECK may never be used directly.
;;;;  Should everything be folded into CALL-WITH-INFO-CHECK?

(define-condition lapack-error (error)
  ;; !! write method for formatting the error message
  ((info :initarg :info :type integer :reader info
	 :documentation "info code")
   (lapack-procedure :initarg :lapack-procedure :type symbol
		     :reader lapack-procedure
		     :documentation "The name without the type prefix (eg 'gesv)."))
  (:documentation "The LAPACK procedure returned a nonzero info
  code."))


(defmacro with-info-check ((procedure-name info-pointer) &body body)
  "Evaluate body with info-pointer bound to an integer, and check it
afterwards, signalling an lapack-error condition if info is nonzero."
  ;;;; !!! how to handle errors nicely? the wrapper can handle this
  ;;;; condition, and output an error message.  There are generally
  ;;;; two kinds of errors in LAPACK: (1) malformed inputs, and (2)
  ;;;; ill-conditioned numerical problems (eg trying to invert a
  ;;;; matrix which is not invertible, etc).
  ;;;;
  ;;;; (1) requires displaying argument names, but in a well-written
  ;;;; library errors like that should not happen, so we will not
  ;;;; worry about it (if it does, that requires debugging & fixing,
  ;;;; not a condition system). In (2), info usually carries all the
  ;;;; information.
  ;;;;
  ;;;; ??? suggestion: an extra argument to the macro on how to
  ;;;; interpret info and generate an error message? create a class
  ;;;; hierarchy of conditions. !!!! do this when API has
  ;;;; stabilized. -- Tamas
  (check-type info-pointer symbol)
  (check-type procedure-name symbol)
  (with-unique-names (info-value)
    `(with-foreign-object (,info-pointer :int32)
       (multiple-value-prog1
	   (progn ,@body)
	 (let ((,info-value (mem-aref ,info-pointer :int32)))
	   (unless (zerop ,info-value)
             ;; ??? two different conditions should be thrown,
             ;; depending on the value of INFO.  Positive INFO usually
             ;; means something substantive (eg a minor is not PSD),
             ;; negative INFO is a bad argument which should never
             ;; happen
             (error 'lapack-error :info ,info-value 
                    :lapack-procedure ',procedure-name)))))))

(defmacro call-with-info-check (procedure &rest arguments)
  "One-liner form that calls the procedure with arguments, and then
checks INFO, which is assumed to be the last argument, using
with-info-check."
  (let ((info-pointer (car (last arguments))))
    (check-type procedure symbol)
    (check-type info-pointer symbol)
    `(with-info-check (,procedure ,info-pointer)
       (funcall ,procedure ,@arguments))))


;;;;  Some (most?) LAPACK procedures allow the caller to query the
;;;;  function for the optimal workspace size, this is a helper macro
;;;;  that does exactly that.  We provide the singular case as a
;;;;  special case if the plural one instead of the other way around,
;;;;  since this would not recurse well.

(defmacro with-work-queries ((&rest specifications) &body body)
  "Call body twice with the given work area specifications, querying
the size for the workspace area.  NOTE: abstraction leaks a bit (body
is there twice), but it should not be a problem in practice.  Body is
most commonly a single function call.

SPECIFICATIONS is a list of triplets (SIZE POINTER LLA-TYPE), where
SIZE and POINTER have to be symbols.  Workspace size (an integer) and
the allocated memory area (pointer) are assigned to these."
  (let* ((sizes (mapcar #'first specifications))
         (pointers (mapcar #'second specifications))
         (returned-sizes (mapcar (lambda (pointer) (gensym* 'returned-size-of- pointer)) pointers))
         (foreign-sizes (mapcar (lambda (pointer) (gensym* 'foreign-size-of- pointer)) pointers))
         (lla-types (mapcar (lambda (pointer) (gensym* 'lla-type-of- pointer)) pointers)))
    (assert (every #'symbolp sizes) () "SIZEs have to be symbols")
    (assert (every #'symbolp pointers) () "POINTERs have to be symbols")
    ;; evaluate lla-types, once only
    `(let ,(mapcar (lambda (lla-type specification)
                     `(,lla-type ,(third specification)))
            lla-types specifications)
       ;; calculate atomic foreign object sizes
       (let ,(mapcar (lambda (foreign-size lla-type)
                       `(,foreign-size (foreign-size* ,lla-type)))
              foreign-sizes lla-types)
         ;; placeholder variables for returned-sizes
         (let ,returned-sizes
           ;; allocate memory for sizes
           (with-foreign-objects ,(mapcar (lambda (size) `(,size :int32 1)) sizes)
             ;; query returned sizes
             ,@(mapcar (lambda (size) `(setf (mem-ref ,size :int32) -1)) sizes)
             (with-foreign-pointers ,(mapcar (lambda (pointer foreign-size)
                                               `(,pointer ,foreign-size))
                                             pointers foreign-sizes)
               ,@body
               ,@(mapcar (lambda (returned-size pointer lla-type size)
                           `(setf ,returned-size (floor (mem-aref* ,pointer ,lla-type))
                                  (mem-ref ,size :int32) ,returned-size))
                         returned-sizes pointers lla-types sizes))
             ;; allocate and call body again
             (with-foreign-pointers ,(mapcar (lambda (pointer foreign-size returned-size)
                                               `(,pointer (* ,foreign-size ,returned-size)))
                                             pointers foreign-sizes returned-sizes)
               ,@body)))))))

(defmacro with-work-query ((size pointer lla-type) &body body)
  "Single-variable version of WITH-WORK-QUERIES."
  `(with-work-queries ((,size ,pointer ,lla-type)) ,@body))
