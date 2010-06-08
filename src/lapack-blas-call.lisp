(in-package :lla)

;;; LAPACK/BLAS functions usually expect either an integer, or one of
;;; the four float LLA types.  The three functions below help
;;; determine the LLA type all the arguments can be coerced to.  It is
;;; always an LLA float type, otherwise an error is signaled.

(defun lb-atom-type% (number)
  "Return the LLA float type that number can be coerced to, or NIL
otherwise.  For (complex) rationals, :(complex)-single is used."
  (typecase number
    (single-float :single)
    (double-float :double)
    (rational :single)
    ((complex single-float) :complex-single)
    ((complex double-float) :complex-double)
    ((complex rational) :complex-single)
    (otherwise (error 'not-within-lla-type))))

(defun lb-vector-type% (vector)
  "Determine LLA type of vector, either by element type, or by
checking each element."
  (aif (array-lla-type vector)
       it
       (reduce #'common-lla-type vector :key #'lb-atom-type%)))

(defun lb-target-type (&rest objects)
  "Find common LLA type objects can be coerced to.  Forces floats
  (from all rational), for use with BLAS/LAPACK.  Not exported."
  (reduce #'common-lla-type objects
          :key (lambda (object)
                 (typecase object
                   (vector (lb-vector-type% object))
                   (elements% (lb-vector-type% (elements object)))
                   (otherwise (lb-atom-type% object))))))

;;; Helper macro that generates the correct LAPACK/BLAS function names
;;; based on a "root" function name.  For some functions (usually
;;; those involving Hermitian matrices), the roots actually differ
;;; based on whether the matrix is real or complex, so the m

(defmacro lb-procedure-name (lla-type name &optional
                              (complex-name name))
  "Evaluate to the LAPACK/BLAS procedure name.  LLA-TYPE has to
evaluate to a symbol denoting a float LLA type.  If you need
conditionals etc, do that outside this macro."
  (check-type name symbol*)
  (check-type complex-name symbol*)
  `(ecase ,lla-type
     (:single #',(make-symbol% "%S" name))
     (:double #',(make-symbol% "%D" name))
     (:complex-single #',(make-symbol% "%C" complex-name))
     (:complex-double #',(make-symbol% "%Z" complex-name))))

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
  ((lapack-procedure :initarg :lapack-procedure :type symbol*
		     :documentation "The name procedure."))
  (:documentation "The LAPACK procedure returned a nonzero info
  code."))

(define-condition lapack-invalid-argument (lapack-error)
  ((position :initarg :position :type fixnum
             :documentation "Position of the illegal argument")
   (argument :initarg :argument :type (or symbol list)
             :documentation "The argument in the call."))
  (:documentation "An argument to a LAPACK procedure had an illegal
  value.  Generally, this indicates a bug in LLA and should not
  happen."))

(define-condition lapack-failure (lapack-error)
  ((index :initarg :index :type fixnum
          :documentation "Index corresponding to error message."))
  (:documentation "Superclass of all LAPACK errors with a positive
  INFO, which is assigned to INDEX."))

(define-condition lapack-singular-matrix (lapack-failure) ())

(defmacro call-with-info-check ((procedure-name info-pointer &optional
                                 condition) &rest arguments)
  "Evaluate body with INFO-POINTER (the last of the arguments by
default) bound to the location of an integer, and check it afterwards.
If INFO is nonzero, signal a LAPACK-INVALID-ARGUMENT for negative
values, and use CONDITION (a subclass of LAPACK-FAILURE, defaults to
the latter) for positive values."
  (unless info-pointer
    (setf info-pointer (car (last arguments))))
  (check-type info-pointer symbol*)
  (check-type procedure-name symbol*)
  (with-unique-names (info-value position)
    (once-only (procedure-name)
      `(with-foreign-object (,info-pointer :int)
         (multiple-value-prog1
             (funcall ,procedure-name ,@arguments)
           (let ((,info-value (mem-aref ,info-pointer :int)))
             (unless (zerop ,info-value)
               (if (minusp ,info-value)
                   (let ((,position (- ,info-value)))
                     (error 'lapack-invalid-argument
                            :lapack-procedure ,procedure-name
                            :argument (nth (1- ,position)
                                           ',arguments)
                            :position ,position))
                   (error ',condition
                          :lapack-procedure ,procedure-name
                          :index ,info-value)))))))))

;;;   We provide the singular case as a
;;; special case if the plural one instead of the other way around,
;;; since this would not recurse well.

(defmacro with-work-queries ((&optional lwork-pointer &rest pointer-type-pairs)
                             &body body)
  "Some (most?) LAPACK procedures allow the caller to query the a function for
the optimal workspace size, this is a helper macro that does exactly that.
LWORK-POINTER is a symbol: a corresponding memory area will be allocated for an
integer and set to -1 before the first pass of BODY.  POINTER-TYPE-PAIRS is a
list of (POINTER TYPE &optional INITIALIZE?), where POINTER is a symbol, and
TYPE evaluates to an LLA type.  In the first pass, a corresponding atom is
allocated, and the result is saved for the second pass, before which a memory
area of the given size is allocated.  The first size is assigned to the memory
are referred to by LWORK-POINTER.  When INITIALIZE?, the atom is initialized
with 0 before the call.  This is useful when the call may not touch the memory
area.

If no arguments are supplied, just return BODY wrapped in a PROGN, once.

NOTE: abstraction leaks a bit (BODY is there twice), but it should not be a
problem in practice.  BODY is most commonly a single function call."
  (unless lwork-pointer
    (return-from with-work-queries `(progn ,@body)))
  (check-type lwork-pointer symbol*)
  (assert (car pointer-type-pairs) () "Need at least one (pointer type) pair.")
  (bind ((pointers (mapcar #'first pointer-type-pairs))
         (type-values (mapcar #'second pointer-type-pairs))
         (initialize? (mapcar #'third pointer-type-pairs))
         ((:flet prepend (prefix names))
          (mapcar (lambda (name) (gensym* prefix name)) names))
         (type-names (prepend '#:type-of- pointers))
         (atom-sizes (prepend '#:atom-size-of- pointers))
         (returned-sizes (prepend '#:returned-size-of- pointers)))
    (assert (every #'symbolp* pointers))
    ;; evaluate types, once only
    `(let ,(mapcar #'list type-names type-values)
       ;; memory size for foreign atoms
       (let ,(mapcar (lambda (atom-size type)
                       `(,atom-size (foreign-size* ,type)))
               atom-sizes type-names)
         ;; placeholders for returned sizes
         (let ,returned-sizes
           ;; allocate memory for lwork
           (with-foreign-object (,lwork-pointer :int)
             ;; first pass - query sizes
             (with-foreign-pointers ,(mapcar #'list pointers atom-sizes)
               ;; set lwork to -1
               (setf (mem-ref ,lwork-pointer :int) -1)
               ;; set first atoms to 0 when needed
               ,@(iter
                   (for pointer :in pointers)
                   (for type-name :in type-names)
                   (for init? :in initialize?)
                   (when init?
                     (collecting `(setf (mem-aref* ,pointer ,type-name)
                                        (zero* ,type-name)))))
               ;; call body
               ,@body
               ;; save sizes
               ,@(mapcar (lambda (returned-size pointer type)
                           `(setf ,returned-size
                                  (as-integer (mem-aref* ,pointer ,type))))
                         returned-sizes pointers type-names)
               ;; also save first one into lwork
               (setf (mem-aref ,lwork-pointer :int) ,(car returned-sizes))
               ;; allocate and call body again
               (with-foreign-pointers
                   ,(mapcar 
                     (lambda (pointer atom-size returned-size)
                       `(,pointer (* ,atom-size ,returned-size)))
                     pointers atom-sizes returned-sizes)
                 ,@body))))))))

;; (defmacro with-work-queries ((&rest specifications) &body body)
;;   "Call body twice with the given work area specifications, querying
;; the size for the workspace area.  NOTE: abstraction leaks a bit (body
;; is there twice), but it should not be a problem in practice.  Body is
;; most commonly a single function call.

;; SPECIFICATIONS is a list of triplets (SIZE POINTER LLA-TYPE), where
;; SIZE and POINTER have to be symbols.  Workspace size (an integer) and
;; the allocated memory area (pointer) are assigned to these."
;;   (unless specifications
;;     (return-from with-work-queries (cons 'progn body)))
;;   (bind ((sizes (mapcar #'first specifications))
;;          (pointers (mapcar #'second specifications))
;;          ((:flet prepend (prefix names))
;;           (mapcar (lambda (name) (gensym* prefix name)) names))
;;          (returned-sizes (prepend '#:returned-size-of- pointers))
;;          (foreign-sizes (prepend '#:foreign-size-of- pointers))
;;          (lla-types (prepend '#:lla-type-of- pointers)))
;;     (assert (every #'symbolp* sizes) () "SIZEs have to be symbols")
;;     (assert (every #'symbolp* pointers) ()
;;             "POINTERs have to be symbols")
;;     ;; evaluate lla-types, once only
;;     `(let ,(mapcar (lambda (lla-type specification)
;;                      `(,lla-type ,(third specification)))
;;             lla-types specifications)
;;        ;; calculate atomic foreign object sizes
;;        (let ,(mapcar (lambda (foreign-size lla-type)
;;                        `(,foreign-size (foreign-size* ,lla-type)))
;;               foreign-sizes lla-types)
;;          ;; placeholder variables for returned-sizes
;;          (let ,returned-sizes
;;            ;; allocate memory for sizes
;;            (with-foreign-objects ,(mapcar (lambda (size)
;;                                             `(,size :int 1)) sizes)
;;              ;; query returned sizes
;;              ,@(mapcar (lambda (size) `(setf (mem-ref ,size :int) -1))
;;                        sizes)
;;              (with-foreign-pointers
;;                  ,(mapcar (lambda (pointer foreign-size)
;;                             `(,pointer ,foreign-size))
;;                           pointers foreign-sizes)
;;                ,@body
;;                ,@(mapcar 
;;                   (lambda (returned-size pointer lla-type size)
;;                     ;; POINTER can be complex, we have to use ABS
;;                     `(setf ,returned-size 
;;                            (as-integer (mem-aref* ,pointer ,lla-type))
;;                            (mem-ref ,size :int) 
;;                            ,returned-size))
;;                   returned-sizes pointers lla-types sizes))
;;              ;; allocate and call body again
;;              (with-foreign-pointers
;;                  ,(mapcar 
;;                    (lambda (pointer foreign-size returned-size)
;;                      `(,pointer (* ,foreign-size ,returned-size)))
;;                    pointers foreign-sizes returned-sizes)
;;                ,@body)))))))

;; (defmacro with-work-query ((size pointer lla-type) &body body)
;;   "Single-variable version of WITH-WORK-QUERIES."
;;   `(with-work-queries% ((,size ,pointer ,lla-type)) ,@body))

(defmacro with-work-area ((pointer lla-type size) &body body)
  "Allocate a work area of size lla-type elements during body,
assigning the pointer to pointer."
  (check-type pointer symbol*)
  `(with-foreign-pointer (,pointer 
			  (* ,size (foreign-size* ,lla-type)))
     ,@body))

(define-with-multiple-bindings with-work-area)


;;;; Miscellaneous utility functions.

(defun zip-eigen% (e-val e-vec check-real? &optional
                   (real-type (array-lla-type e-val)))
  "Some LAPACK routines (yes, xGEEV, I am talking to you) return real and
complex parts of eigenvalues in two separate vectors.  LLA routines will try to
collect these in a single vector that is the concatenation of vectors containing
the real and imaginary parts, and this function will make a single vector out of
them, given as E-VAL.  If E-VEC is non-nil, it is taken to be a corresponding
eigenvector matrix, and eigenvectors are collected from it.  If CHECK-REAL?,
eigenvectors are returned as real numbers (instead of complex) when all of them
are real, otherwise they are always complex.  The second value returned is the
matrix of eigenvectors."
  ;; !!! could do with more optimization and benchmarking
  (declare (optimize speed (safety 0)))
  (macrolet ((zip-eigen%% (real-type)
               (let ((complex-type (complex-lla-type real-type)))
                 `(let* ((n (aif (divides? (length e-vec) 2) it
                                 (error "E-VEC has to contain an even number ~
                                         of elements")))
                         (val-zipped (lla-vector ,complex-type n))
                         (vec-zipped (when e-vec
                                       (let ((n^2 (expt n 2)))
                                         (assert (= (length e-vec) n^2))
                                         (lla-vector ,complex-type n^2)))))
                    (declare (fixnum n)
                             (type ,(lla-vector-type real-type) e-val e-vec)
                             (type ,(lla-vector-type complex-type) val-zipped))
                    (iter
                      (with index := 0)
                      (with all-real? := t)
                      (declare (iterate:declare-variables))
                      (let* ((realpart (aref e-val index))
                             (imagpart (aref e-val (+ index n)))
                             (real? (zerop imagpart)))
                        (if real?
                            (progn
                              ;; copy eigenvalue
                              (setf (aref val-zipped index)
                                    (complex realpart 0))
                              ;; copy a single real column
                              (when e-vec
                                (copy-elements e-vec (cm-index2 n 0 index)
                                               vec-zipped (cm-index2 n 0 index) n))
                              ;; increment index
                              (incf index))
                            (progn
                              ;; copy eigenvalue and conjugate
                              (setf (aref val-zipped index)
                                    (complex realpart imagpart)
                                    (aref val-zipped (1+ index))
                                    (complex realpart (- imagpart)))
                              ;; copy two compex columns, conjugates of each other
                              (let* ((vec-index-real (cm-index2 n 0 index))
                                     (vec-index-imag (+ vec-index-real n)))
;                                (declare (fixnum vec-index-real vec-index-imag))
                                (iter
                                  (for (the fixnum real-index)
                                       :from vec-index-real
                                       :to vec-index-imag)
                                  (for (the fixnum imag-index) :from vec-index-imag)
                                  (declare (iterate:declare-variables))
                                  (let ((realpart (aref e-vec real-index))
                                        (imagpart (aref e-vec imag-index)))
                                    (setf (aref vec-zipped real-index)
                                          (complex realpart imagpart)
                                          (aref vec-zipped imag-index)
                                          (complex realpart (- imagpart))))))
                              ;; increment index
                              (incf index 2)))
                        ;; termination
                        (while (< index n))
                        (when (and all-real? (not real?))
                          (setf all-real? nil))
                        (finally
                         (flet ((conv (vector)
                                  (if (and check-real? all-real?)
                                      (map (lla-vector-type real-type) #'realpart vector)
                                      vector)))
                           (return (values (conv val-zipped)
                                           (when vec-zipped
                                             (make-matrix% n n 
                                                           (conv vec-zipped)))))))))))))
    (ecase real-type
      (:single (zip-eigen%% :single))
      (:double (zip-eigen%% :double)))))

(defun zip-eigenvalues (val-pointer n real-type complex-type
                        check-real-p)
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
                       (always (zerop (mem-aref* val-pointer 
                                                 real-type 
                                                 (+ n i))))))))
    (if real-p
        ;; no complex eigenvalues
        (let ((elements (lla-vector real-type n)))
          (iter
            (for i :from 0 :below n)
            (setf (aref elements i) 
                  (mem-aref* val-pointer real-type i)))
          (values elements nil))
        ;; complex eigenvalues
        (let ((elements (lla-vector complex-type n)))
          (iter
            (for i :from 0 :below n)
            (setf (aref elements i) 
                  (complex (mem-aref* val-pointer real-type i)
                           (mem-aref* val-pointer real-type (+ n i)))))
          (values elements t)))))

(defun zip-eigenvectors (val-pointer vec-pointer n
                         real-type complex-type)
  "Collect complex eigenvectors from S/DGEEV.  Should only be called
when the second value returned by ZIP-EIGENVALUES is non-nil."
  (prog ((vec (lla-vector complex-type (* n n)))
         (i 0))
   top
     (let ((column-start-index (cm-index2 n 0 i)))
       (if (zerop (mem-aref* val-pointer real-type (+ n i)))
           ;; real
           (progn
             (iter
               (for j :from 0 :below n)
               (for vec-index :from column-start-index)
               (setf (aref vec vec-index)
                     (complex 
                      (mem-aref* vec-pointer real-type vec-index))))
             (incf i))
           ;; complex, assemble from real +- imaginary columns
           (progn
             (iter
               (for j :from 0 :below n)
               (for vec-index :from column-start-index)
               (for vec-index2 :from (+ column-start-index n))
               (with realpart := 
                     (mem-aref* vec-pointer real-type vec-index))
               (with imagpart := 
                     (mem-aref* vec-pointer real-type vec-index2))
               (setf (aref vec vec-index) (complex realpart imagpart)
                     (aref vec vec-index2) (complex realpart 
                                                    (- imagpart))))
             (incf i 2))))
     (if (< i n)
         (go top)
         (return vec))))

;;; Collecting the matrix/vector at the end.

(defun sum-last-rows (elements m nrhs n)
  "Sum & return (as a NUMERIC-VECTOR) the last (- M N) rows of an M x
NRHS matrix, given as a Lisp vector in column-major view.  NOTE:
needed to interface to LAPACK routines like xGELS."
  (bind ((row-sums (lla-vector (real-lla-type 
                                (array-lla-type elements))
                               nrhs)))
    (dotimes (col nrhs)
      (setf (aref row-sums col)
            (sse-elements% elements
                           (cm-index2 m n col)
                           (cm-index2 m m col))))
    row-sums))

;;;; nice interface for matrices, probably the most important macro


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

;;; lb-call, the main interface macro

(defclass lb-bindings ()
  ((bindings :initform nil :documentation "Filtered bindings list.")
   (atoms :initform nil)
   (vectors :initform nil)
   (set-restricted? :initform nil :documentation "List of 
  (name set-restricted?) pairs on which to call SET-RESTRICTED.")
   (outputs :initform nil)
   (info :initform nil :documentation
         "When set, (info-pointer condition)")
   (work-queries :initform nil
                 :documentation "Will be passed to WITH-WORK-QUERIES.")
   (work-areas :initform nil)))

(defun lb-process-binding (binding-form bindings)
  "Process a binding form, putting the result in bindings."
  (bind (((variable-form &rest value-forms) binding-form))
    (if (atom variable-form)
        (push `(,variable-form ,@value-forms)
              (slot-value bindings 'bindings))
        (lb-process-binding% (car variable-form)
                             (cdr variable-form)
                             value-forms
                             bindings))))

(defmacro lb-collect ((bindings) &body slot-value-pairs)
  "Helper function for use inside lb-process-binding% methods.

   Example:

     (lb-collect (bindings)
       bindings `(,name ,value)
       atoms (make-instance 'lb-binding-atom ...)
       process `((:atom ...) value))"
  (assert (divides? (length slot-value-pairs) 2))
  (once-only (bindings)
    `(progn
       ,@(iter
           (for (slot value &rest rest) :on slot-value-pairs
                :by #'cddr)
           (collecting 
             (if (eq slot 'process)
                 `(lb-process-binding ,value ,bindings)
                 `(push ,value (slot-value ,bindings ',slot))))))))

(defgeneric lb-process-binding% (keyword specification value-forms
                                 bindings)
  (:documentation "Return filtered bindings, add forms to BINDINGS."))

(defmethod lb-process-binding% (keyword specification value-forms
                                bindings)
  ;; will be processed as is
  (lb-collect (bindings)
    bindings `((,keyword ,@specification) ,@value-forms)))

(defmethod lb-process-binding% ((keyword (eql :atom))
                                specification value-forms bindings)
  ;; syntax: ((:atom pointer-name lla-type) value)
  (bind (((pointer-name lla-type-value) specification)
         (lla-type-name (gensym* '#:lla-type- pointer-name))
         (value-name (gensym* '#:atom- pointer-name)))
    (lb-collect (bindings)
      bindings `(,value-name ,@value-forms)
      bindings `(,lla-type-name ,lla-type-value)
      atoms `(,pointer-name ,lla-type-name ,value-name))))

(defmethod lb-process-binding% ((keyword (eql :char))
                                specification value-forms bindings)
  ;; syntax: ((:char pointer-name) value), expands to
  ;; (:atom pointer-name :char)
  (bind (((pointer-name) specification))
    (check-type pointer-name symbol*)
    (lb-collect (bindings)
      process `((:atom ,pointer-name :char) ,@value-forms))))

(defmethod lb-process-binding% ((keyword (eql :integer))
                                specification value-forms bindings)
  ;; syntax: ((:integer pointer-name) value), expands to 
  ;; (:atom pointer-name :integer)
  (bind (((pointer-name) specification))
    (check-type pointer-name symbol*)
    (lb-collect (bindings)
      process `((:atom ,pointer-name :integer) ,@value-forms))))

(defmethod lb-process-binding% ((keyword (eql :vector))
                                specification value-forms bindings)
  ;; syntax: ((:vector pointer-name lla-type &optional copy-or-output) vector)
  (bind (((pointer-name lla-type-value &optional
                        output) specification)
         (lla-type-name (gensym* '#:lla-type- pointer-name))
         (value-name (gensym* '#:vector- pointer-name)))
    (check-type pointer-name symbol*)
    (check-type output (or null (eql :copy) symbol*))
    (lb-collect (bindings)
      bindings `(,value-name ,@value-forms)
      bindings `(,lla-type-name ,lla-type-value)
      vectors `(,value-name ,pointer-name ,lla-type-name ,output))
    (awhen (and (not (eq output :copy)) output)
      (check-type it symbol*)
      (lb-collect (bindings) bindings it))))

(defmethod lb-process-binding% ((keyword (eql :output))
                                specification value-forms bindings)
  ;; syntax: ((:output pointer-name lla-type output) length)
  (bind (((pointer-name lla-type-value output-name)
          specification)
         (lla-type-name (gensym* '#:lla-type- pointer-name))
         (length-name (gensym* '#:length- pointer-name)))
    (check-type pointer-name symbol*)
    (check-type output-name symbol*)
    (lb-collect (bindings)
      bindings `(,length-name ,@value-forms)
      bindings `(,lla-type-name ,lla-type-value)
      outputs `(,output-name ,pointer-name ,lla-type-name ,length-name)
      bindings output-name)))

(defmethod lb-process-binding% ((keyword (eql :dimension))
                                specification value-forms bindings)  
  ;; syntax: ((:dimension &optional var pointer) value)
  ;; dim-spec: NIL (not processed), VAR, (VAR), or (VAR POINTER).
  (when specification
    (bind (((variable-name &optional pointer-name) specification))
      (lb-collect (bindings)
        bindings `(,variable-name ,@value-forms))
      (when pointer-name
        (lb-collect (bindings)
          atoms `(,pointer-name :integer ,variable-name))))))

(defmethod lb-process-binding% ((keyword (eql :matrix))
                                specification value-forms bindings)
  ;; syntax: ((:matrix pointer lla-type nrow-spec ncol-spec &key
  ;;             (set-restricted? t) output) matrix)
  ;; nrow/ncol-spec: see (:dimension ...)
  (bind (((pointer-name lla-type-value nrow-spec ncol-spec
                        &key (set-restricted? t)
                        output) specification)
         (value-name (gensym* '#:value- pointer-name)))
    (lb-collect (bindings)
      bindings `(,value-name ,@value-forms)
      process `((:vector ,pointer-name ,lla-type-value
                         ,output) (elements ,value-name))
      process `((:dimension ,@(mklist nrow-spec))
                (nrow ,value-name))
      process `((:dimension ,@(mklist ncol-spec))
                (ncol ,value-name)))
    (when set-restricted?
      (lb-collect (bindings)
        set-restricted? `(,value-name ,set-restricted?)))))

(defmethod lb-process-binding% ((keyword (eql :check))
                                specification value-forms bindings)
  ;; syntax: ((:check info-pointer &optional (condition lapack-failure)))
  (bind (((:slots info) bindings)
         ((info-pointer &optional (condition 'lapack-failure))
          specification))
    (assert (not info) () ":CHECK already specified.")
    (check-type info-pointer symbol*)
    (check-type value-forms null)
    (lb-collect (bindings)
      info `(,info-pointer ,condition))))

(defmethod lb-process-binding% ((keyword (eql :work-queries))
                                specification value-forms bindings)
  ;; syntax: ((:work-queries lwork% (work% type) ...)), where lwork% is the
  ;; pointer to the size of the work area, work% to the work area, and type is
  ;; their LLA type. !!! specifications is not expanded, just don't use any side
  ;; effects or depend on the order.
  (bind (((:slots work-queries) bindings))
    (assert (not work-queries) () ":WORK-QUERIES already specified.")
    (assert specification () ":WORK-QUERIES needs a non-empty specification.")
    (assert (not value-forms) () ":WORK-QUERIES does not take value forms.")
    (setf work-queries specification)))

(defmethod lb-process-binding% ((keyword (eql :work))
                                specification value-forms bindings)
  ;; syntax: ((:work pointer lla-type) size)
  (bind (((pointer-name lla-type-value) specification)
         (lla-type-name (gensym* '#:lla-type- pointer-name))
         (size-name (gensym* '#:lla-type- pointer-name)))
    (check-type pointer-name symbol*)
    (lb-collect (bindings)
      bindings `(,lla-type-name ,lla-type-value)
      bindings `(,size-name ,@value-forms)
      work-areas `(,pointer-name ,lla-type-name ,size-name))))


(defmacro lb-call (binding-forms &body body)
  "The BINDING-FORMs below are captured, the rest are passed to BIND
as is.  See the comments in the corresponding LB-PROCESS-BINDING%
method.

 ((:atom pointer lla-type) value)

 ((:char pointer) value) [shorthand for ((:atom pointer :char) ...)]

 ((:vector pointer lla-type &optional copy-or-output) vector)

 ((:output pointer-name lla-type output) length)

 ((:matrix pointer lla-type nrow-spec ncol-spec &key
     (set-restricted? t) output) matrix)

 ((:check info-pointer &optional (condition lapack-failure)))

 ((:work-queries lwork% (work% type) ...))

 ((:work pointer type) size)

The following are meant for internal use:

 ((:dimension &optional var pointer) value)

"
  (bind ((bindings (make-instance 'lb-bindings))
         ((:flet expand (slot-name))
          (nreverse (slot-value bindings slot-name))))
    (dolist (binding-form binding-forms)
      (lb-process-binding binding-form bindings))
    `(bind ,(nreverse (slot-value bindings 'bindings))
       ;; set restricted elements
       ,@(mapcar (lambda (sr)
                   `(when ,(second sr)
                      (set-restricted ,(first sr))))
                 (expand 'set-restricted?))
       ;; define CALL macro
       (macrolet ((call (procedure &rest arguments)
                    `(with-fortran-atoms ,',(expand 'atoms)
                       (with-pinned-vectors ,',(expand 'vectors)
                         (with-vector-outputs ,',(expand 'outputs)
                           (with-work-queries ,',(slot-value bindings
                                                             'work-queries)
                             (with-work-areas ,',(expand 'work-areas)
                               ,,(bind (((:slots info) bindings))
                                   (if info
                                       ``(call-with-info-check (,procedure
                                                                ,@',@info)
                                                               ,@arguments)
                                       ``(funcall ,procedure ,@arguments))))))))))
         ,@body))))
