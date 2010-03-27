(in-package :lla)

;;;; Code for interfacing with Fortran calls.


(declaim (inline lb-target-type))
(defun lb-target-type (&rest objects)
  "Find common target type of objects to.  Forces floats, should be
used in LAPACK."
  (binary-code->lla-type
   (reduce #'logior objects :key (lambda (object)
                                   (lla-type->binary-code
                                    (lla-type object))))))

;;; Helper functions that generate the correct LAPACK/BLAS function
;;; names based on a "root" function name.  For some functions
;;; (usually those involving Hermitian matrices), the roots actually
;;; differ based on whether the matrix is real or complex.

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

;;; Some LAPACK procedures can signal errors, they do this via an INFO
;;; output integer.  Here we provide macros to capture these errors.
;;;
;;; Currently, a general lapack-error condition is thrown.  This may
;;; be too coarse, as INFO actually provides quite a bit of
;;; information of what went wrong.  There are usually two kinds of
;;; errors: invalid/nonsensical arguments (should never happen), and
;;; incomputable problems (matrix being singular, usually with
;;; detailed info).  !!! We should capture the latter and provide more
;;; sensible error messages.  Maybe instead of throwing LAPACK-ERROR,
;;; macros should accept a form that tells what kind of error to
;;; signal.  Expect CALL-WITH-INFO-CHECK to be modified in the future,
;;; eg all current arguments go into a list, then body gets PROCEDURE
;;; and INFO values in case of an error and may do what it wants with
;;; them.
;;;
;;; ??? Most (maybe all?) LAPACK functions just have INFO as their
;;; last argument, so WITH-INFO-CHECK may never be used directly.
;;; Should everything be folded into CALL-WITH-INFO-CHECK?

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
  ;; !!! how to handle errors nicely? the wrapper can handle this
  ;; condition, and output an error message.  There are generally
  ;; two kinds of errors in LAPACK: (1) malformed inputs, and (2)
  ;; ill-conditioned numerical problems (eg trying to invert a
  ;; matrix which is not invertible, etc).
  ;;
  ;; (1) requires displaying argument names, but in a well-written
  ;; library errors like that should not happen, so we will not
  ;; worry about it (if it does, that requires debugging & fixing,
  ;; not a condition system). In (2), info usually carries all the
  ;; information.
  ;;
  ;; ??? suggestion: an extra argument to the macro on how to
  ;; interpret info and generate an error message? create a class
  ;; hierarchy of conditions. !!!! do this when API has
  ;; stabilized. -- Tamas
  (check-type info-pointer symbol)
  (check-type procedure-name symbol)
  (with-unique-names (info-value)
    `(with-foreign-object (,info-pointer :int)
       (multiple-value-prog1
	   (progn ,@body)
	 (let ((,info-value (mem-aref ,info-pointer :int)))
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


;;; Some (most?) LAPACK procedures allow the caller to query the
;;; function for the optimal workspace size, this is a helper macro
;;; that does exactly that.  We provide the singular case as a
;;; special case if the plural one instead of the other way around,
;;; since this would not recurse well.

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
           (with-foreign-objects ,(mapcar (lambda (size) `(,size :int 1)) sizes)
             ;; query returned sizes
             ,@(mapcar (lambda (size) `(setf (mem-ref ,size :int) -1)) sizes)
             (with-foreign-pointers ,(mapcar (lambda (pointer foreign-size)
                                               `(,pointer ,foreign-size))
                                             pointers foreign-sizes)
               ,@body
               ,@(mapcar (lambda (returned-size pointer lla-type size)
                           ;; POINTER can be complex, we have to use ABS too
                           `(setf ,returned-size (floor (abs (mem-aref* ,pointer ,lla-type)))
                                  (mem-ref ,size :int) ,returned-size))
                         returned-sizes pointers lla-types sizes))
             ;; allocate and call body again
             (with-foreign-pointers ,(mapcar (lambda (pointer foreign-size returned-size)
                                               `(,pointer (* ,foreign-size ,returned-size)))
                                             pointers foreign-sizes returned-sizes)
               ,@body)))))))

(defmacro with-work-query ((size pointer lla-type) &body body)
  "Single-variable version of WITH-WORK-QUERIES."
  `(with-work-queries ((,size ,pointer ,lla-type)) ,@body))

(defmacro with-work-area ((pointer lla-type size) &body body)
  "Allocate a work area of size lla-type elements during body,
assigning the pointer to pointer."
  (check-type pointer symbol)
  `(with-foreign-pointer (,pointer 
			  (* ,size (foreign-size* ,lla-type)))
     ,@body))

(define-with-multiple-bindings with-work-area)


;;;; Miscellaneous utility functions.

(defun zip-eigenvalues (val-pointer n real-type complex-type check-real-p)
  "Return the complex numbers stored at VAL-POINTER (N real parts,
followed by N imaginary parts) as a NUMERIC-VECTOR (either
SINGLE/DOUBLE or COMPLEX-SINGLE/COMPLEX-DOUBLE).  If CHECK-REAL-P,
then check if the imaginary part is 0 and if so, return a
NUMERIC-VECTOR-SINGLE/DOUBLE, otherwise always return a complex one.
The second value is non-nil if there are complex eigenvalues.  *Usage
note:* some LAPACK routines return real and imaginary parts of vectors
separately, we have to assemble them. *NOT EXPORTED*."
  (let ((real-p (and check-real-p 
                     (iter
                       (for i :from 0 :below n)
                       (always (zerop (mem-aref val-pointer 
                                                real-type (+ n i))))))))
    (if real-p
        ;; no complex eigenvalues
        (let ((elements (make-nv-elements real-type n)))
          (iter
            (for i :from 0 :below n)
            (setf (aref elements i) (mem-aref val-pointer real-type i)))
          (values (make-nv* real-type elements) nil))
        ;; complex eigenvalues
        (let ((elements (make-nv-elements complex-type n)))
          (iter
            (for i :from 0 :below n)
            (setf (aref elements i) (complex (mem-aref val-pointer
                                                       real-type i)
                                             (mem-aref val-pointer
                                                       real-type (+ n i)))))
          (values (make-nv* complex-type elements) t)))))

(defun zip-eigenvectors (val-pointer vec-pointer n real-type complex-type)
  "Collect complex eigenvectors from S/DGEEV.  Should only be called
when the second value returned by ZIP-EIGENVALUES is non-nil."
  (prog ((vec (make-nv-elements complex-type (* n n)))
         (i 0))
   top
     (let ((column-start-index (cm-index2 n 0 i)))
       (if (zerop (mem-aref val-pointer real-type (+ n i)))
           ;; real
           (progn
             (iter
               (for j :from 0 :below n)
               (for vec-index :from column-start-index)
               (setf (aref vec vec-index)
                     (complex (mem-aref vec-pointer real-type vec-index))))
             (incf i))
           ;; complex, assemble from real +- imaginary columns
           (progn
             (iter
               (for j :from 0 :below n)
               (for vec-index :from column-start-index)
               (for vec-index2 :from (+ column-start-index n))
               (with realpart := (mem-aref vec-pointer real-type vec-index))
               (with imagpart := (mem-aref vec-pointer real-type vec-index2))
               (setf (aref vec vec-index) (complex realpart imagpart)
                     (aref vec vec-index2) (complex realpart (- imagpart))))
             (incf i 2))))
     (if (< i n)
         (go top)
         (return vec))))

;;; Collecting the matrix/vector at the end.

(defmacro ifor ((variable from to &optional (name nil))
                &body body)
  "FOR loop for integers.  Declares the variables to be fixnums."
  (check-type variable symbol)
  (check-type name symbol)
  (let ((loop-top-name (make-symbol* 'loop-top- name))
        (end-index (gensym "END")))
    `(let ((,variable ,from)
           (,end-index ,to))
       (declare (fixnum ,variable ,end-index))
       (block ,name
         (tagbody
         ,loop-top-name
            ,@body
            (incf ,variable)
            (when (< ,variable ,end-index)
              (go ,loop-top-name)))))))

(defgeneric sum-last-rows (vector lla-type m nrhs n)
  (:documentation "Sum & return (as a NUMERIC-VECTOR) the last m-n
rows of an m x nrhs matrix, given as a Lisp vector in column-major
view.  NOTE: needed to interface to LAPACK routines like xGELS."))
(expand-for-lla-types (lla-type :exclude-integer-p t)
  (bind (((result-type complex-p) (ecase lla-type
                                    (:single '(:single nil))
                                    (:double '(:double nil))
                                    (:complex-single '(:single t))
                                    (:complex-double '(:double t)))))
    `(defmethod sum-last-rows (vector (lla-type (eql ,lla-type)) m nrhs n)
       (declare (optimize (debug 3)) ; (optimize speed (safety 0))
                (type ,(nv-array-type lla-type) vector)
                (type fixnum m nrhs n))
       (let* ((result (make-nv-elements ,result-type nrhs))
              (nrow (- m n)))
         (dotimes (col nrhs)
           (declare (fixnum col))
           (let ((sum ,(coerce 0 (lla-type->lisp-type result-type)))
                 (start-index (the fixnum (+ (the fixnum (* col m)) n))))
             (ifor (vector-index start-index (the fixnum (+ start-index nrow)))
               (incf sum ,(if complex-p
                              ;; here we need the square *modulus*
                              `(let ((x (aref vector vector-index)))
                                 (+ (expt (realpart x) 2) (expt (imagpart x) 2)))
                              `(expt (aref vector vector-index) 2))))
             (setf (aref result col) sum)))
         (make-nv* ,result-type result)))))

;;;; nice interface for matrices, probably the most important macro

(defmacro with-matrix-input (((matrix nrow ncol leading-dimension &optional
                                      keyword output)
                              pointer lla-type
                              &optional (set-restricted t)) &body body)
  "Convenience macro for using matrices as input for Fortran calls.
Essentially like WITH-NV-INPUT, except that matrix dimensions are
bound to variables.  NROW and NCOL have the following syntax:
  variable-name  -- dimensions is assigned to the variable
  (variable-name fortran-atom-name)  --  value is also assigned to a
     memory location as an integer using WITH-FORTRAN-ATOM."
  (flet ((parse-dimspec (dimspec)
           "Parse dimension name binding specification."
           (cond
             ((and dimspec (atom dimspec))
              (check-type dimspec symbol)
              (list dimspec nil))
             ((and dimspec (listp dimspec) (= (length dimspec) 2))
              (bind (((name fortran-atom-name) dimspec))
                (assert (every #'symbolp dimspec) ()
                        "Need symbol for dimension specification.")
                `(,name ((:integer ,fortran-atom-name ,name)))))
             (t (error "Invalid dimension specification ~A.  ~
                       Use NAME or (NAME FORTRAN-ATOM-NAME)." dimspec)))))
    (bind (((nrow-name nrow-fortran-atom-expansion) (parse-dimspec nrow))
           ((ncol-name ncol-fortran-atom-expansion) (parse-dimspec ncol)))
      (once-only (matrix)
        `(bind (((:slots-read-only (,nrow-name nrow) (,ncol-name ncol)) ,matrix))
           (if ,set-restricted
               (set-restricted ,matrix))
           (with-fortran-atoms (,@nrow-fortran-atom-expansion
                                ,@ncol-fortran-atom-expansion
                                (:integer ,leading-dimension (leading-dimension ,matrix)))
             (with-nv-input ((,matrix ,keyword ,output) ,pointer ,lla-type)
               ,@body)))))))

(define-with-multiple-bindings with-matrix-input)

;;; floating point traps
;;;
;;; Apparently, the only trap that we need to mask is division by
;;; zero, and that only for a few operations.  Non-numerical floating
;;; points values are used internally (eg in SVD calculations), but
;;; only reals are returned.

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
