(in-package :lla)

(defun lb-transpose (transpose library &optional complex?)
  "Return an integer (enum) that denotes the desired transpose operation (NIL:
no transpose, T: transpose, *: conjugate transpose which is usually equivalent
to T for real-valued matrices.  When complex?, T is replaced by *."
  (when (and transpose complex?)
    (setf transpose '*))
  (ecase library
    (:lapack
       (ecase transpose
         ((nil) +n+)
         ((t) +t+)
         ((*) +c+)))
    (:blas
       (ecase transpose
         ((nil) :CblasNoTrans)
         ((t) :CblasTrans)
         ((*) :CblasConjTrans)))))

(defun transposed-dimensions (m n transpose &optional (library :blas))
  "Process dimensions of a matrix, which may be transposed.  Third value is
the enum expected by the library."
  (ecase transpose
    ((nil) (values m n (ecase library (:blas ) (:lapack (char-code #\N)))))
    ((t) (values n m (ecase library (:blas) (:lapack #\T))))
    ((*) (values n m (ecase library (:blas ) (:lapack #\C))))))

(defun maybe-vector-as-matrix (vector-or-matrix orientation &optional 
                               transpose (library :blas))
  "Process dimensions of vectors which are meant to be reoriented as matrices.
Orientation is :ROW or :COLUMN, and TRANSPOSE indicates transposition (see
TRANS).  Return (values DIMENSION0 DIMENSION1 ORIENTATION LEADING-DIMENSION
TRANSPOSE-ENUM).  ORIENTATION is NIL or the ORIENTATION argument to the
function, depending on whether the argument was a vector."
  (let+ (((d1 &optional d2 d-rest) (array-dimensions vector-or-matrix))
         ((&flet error% () (error "~A is not a vector or a matrix."
                                  vector-or-matrix)))
         ((&values m n orientation)
          (cond
            (d-rest (error%))
            (d2 (values d1 d2))
            (d1 (ecase orientation
                  (:row (values 1 d1 :row))
                  (:column (values d1 1 :column))))
            (t (error%))))
         ((&values m-trans n-trans) (if transpose
                                        (values n m)
                                        (values m n))))
    (values m-trans n-trans orientation n (lb-transpose transpose library))))

(defun vector-or-matrix-dimensions (m n orientation)
  "Return dimensions for matrices that can potentially be vectors, depending
 on orientation."
  (ecase orientation
    (:row (assert (= m 1)) n)
    (:column (assert (= n 1)) m)
    ((nil) (list m n))))

(defun maybe-pick-first-element (array pick?)
  "When PICK?, return first element of array, otherwise the whole array."
  (if pick?
      (progn (assert (= (array-total-size array) 1) ()
                     "~A is supposed to have only one element." array)
             (row-major-aref array 0))
      array))

(defun matrix-from-first-rows (matrix nrow orientation)
  "Create a matrix (or vector, depending on ORIENTATION) from the first rows
NRHS of MATRIX.  Used for interfacing with xGELS, extracting R from QR
decompositions, etc."
  (let+ (((nil n) (array-dimensions matrix)))
    (copy-array (displace-array matrix 
                                (vector-or-matrix-dimensions nrow n
                                                             orientation)))))
(defparameter *lla-double?* t
  "Determines whether rational->float conversions result in double or single
   floats.")

(defun common-float-type (&rest objects)
  "Determine common float type for OBJECTS.  For use in LAPACK/BLAS calls."
  (common-lla-type objects :force-float? t :double? *lla-double?*))

(defmacro procedure-name (library lla-type name &optional
                             (complex-name name))
  "Evaluate to the LAPACK/BLAS procedure name.  LIBRARY is :BLAS or :LAPACK.
LLA-TYPE has to evaluate to a symbol denoting a float LLA type.  If you need
conditionals etc for the function name, do that outside this macro.  For some
functions (usually those involving Hermitian matrices), the names actually
differ based on whether the matrix is real or complex, use COMPLEX-NAME in
that case."
  (check-type name symbol*)
  (check-type complex-name symbol*)
  (let ((library-prefix (ecase library
                          (:blas '#:CBLAS_)
                          (:lapack '#:LAPACKE_))))
    `(ecase ,lla-type
       (:single #',(make-symbol% library-prefix "S" name))
       (:double #',(make-symbol% library-prefix "D" name))
       (:complex-single #',(make-symbol% library-prefix "C" complex-name))
       (:complex-double #',(make-symbol% library-prefix "Z" complex-name)))))

(defun matrix-layout (library layout)
  "Return the matrix layout constant."
  (ecase library
    (:lapack (ecase layout
               (:row-major LAPACK_ROW_MAJOR)
               (:column-major LAPACK_COL_MAJOR)
               ((nil) nil)))
    (:blas (ecase layout
               (:row-major :CBLASROWMAJOR)
               (:column-major :CBLASCOLMAJOR)
               ((nil) nil)))))

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
		     :documentation "The name of the procedure."))
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
  ((info :initarg :info :type fixnum
          :documentation "INFO corresponding to error message."))
  (:documentation "Superclass of all LAPACK errors with a positive INFO"))

(define-condition lapack-singular-matrix (lapack-failure) ())

;; (defmacro call-with-info-check ((procedure-name info-pointer &optional
;;                                  condition) &rest arguments)
;;   "Evaluate body with INFO-POINTER (the last of the arguments by
;; default) bound to the location of an integer, and check it afterwards.
;; If INFO is nonzero, signal a LAPACK-INVALID-ARGUMENT for negative
;; values, and use CONDITION (a subclass of LAPACK-FAILURE, defaults to
;; the latter) for positive values."
;;   (unless info-pointer
;;     (setf info-pointer (car (last arguments))))
;;   (check-type info-pointer symbol*)
;;   (check-type procedure-name symbol*)
;;   (with-unique-names (info-value position)
;;     (once-only (procedure-name)
;;       `(with-foreign-object (,info-pointer :int)
;;          (multiple-value-prog1
;;              (funcall ,procedure-name ,@arguments)
;;            (let ((,info-value (mem-aref ,info-pointer :int)))
;;              (unless (zerop ,info-value)
;;                (if (minusp ,info-value)
;;                    (let ((,position (- ,info-value)))
;;                      (error 'lapack-invalid-argument
;;                             :lapack-procedure ,procedure-name
;;                             :argument (nth (1- ,position)
;;                                            ',arguments)
;;                             :position ,position))
;;                    (error ',condition
;;                           :lapack-procedure ,procedure-name
;;                           :index ,info-value)))))))))

;; (defmacro with-work-queries ((&optional lwork-pointer &rest pointer-type-pairs)
;;                              &body body)
;;   "Some (most?) LAPACK procedures allow the caller to query the a function for
;; the optimal workspace size, this is a helper macro that does exactly that.
;; LWORK-POINTER is a symbol: a corresponding memory area will be allocated for an
;; integer and set to -1 before the first pass of BODY.  POINTER-TYPE-PAIRS is a
;; list of (POINTER TYPE &optional INITIALIZE?), where POINTER is a symbol, and
;; TYPE evaluates to an LLA type.  In the first pass, a corresponding atom is
;; allocated, and the result is saved for the second pass, before which a memory
;; area of the given size is allocated.  The first size is assigned to the memory
;; are referred to by LWORK-POINTER.  When INITIALIZE?, the atom is initialized
;; with 0 before the call.  This is useful when the call may not touch the memory
;; area.

;; If no arguments are supplied, just return BODY wrapped in a PROGN, once.

;; NOTE: abstraction leaks a bit (BODY is there twice), but it should not be a
;; problem in practice.  BODY is most commonly a single function call."
;;   (unless lwork-pointer
;;     (return-from with-work-queries `(progn ,@body)))
;;   (check-type lwork-pointer symbol*)
;;   (assert (car pointer-type-pairs) () "Need at least one (pointer type) pair.")
;;   (bind ((pointers (mapcar #'first pointer-type-pairs))
;;          (type-values (mapcar #'second pointer-type-pairs))
;;          (initialize? (mapcar #'third pointer-type-pairs))
;;          ((:flet prepend (prefix names))
;;           (mapcar (lambda (name) (gensym* prefix name)) names))
;;          (type-names (prepend '#:type-of- pointers))
;;          (atom-sizes (prepend '#:atom-size-of- pointers))
;;          (returned-sizes (prepend '#:returned-size-of- pointers)))
;;     (assert (every #'symbolp* pointers))
;;     ;; evaluate types, once only
;;     `(let ,(mapcar #'list type-names type-values)
;;        ;; memory size for foreign atoms
;;        (let ,(mapcar (lambda (atom-size type)
;;                        `(,atom-size (foreign-size* ,type)))
;;                atom-sizes type-names)
;;          ;; placeholders for returned sizes
;;          (let ,returned-sizes
;;            ;; allocate memory for lwork
;;            (with-foreign-object (,lwork-pointer :int)
;;              ;; first pass - query sizes
;;              (with-foreign-pointers ,(mapcar #'list pointers atom-sizes)
;;                ;; set lwork to -1
;;                (setf (mem-ref ,lwork-pointer :int) -1)
;;                ;; set first atoms to 0 when needed
;;                ,@(iter
;;                    (for pointer :in pointers)
;;                    (for type-name :in type-names)
;;                    (for init? :in initialize?)
;;                    (when init?
;;                      (collecting `(setf (mem-aref* ,pointer ,type-name)
;;                                         (zero* ,type-name)))))
;;                ;; call body
;;                ,@body
;;                ;; save sizes
;;                ,@(mapcar (lambda (returned-size pointer type)
;;                            `(setf ,returned-size
;;                                   (as-integer (mem-aref* ,pointer ,type))))
;;                          returned-sizes pointers type-names)
;;                ;; also save first one into lwork
;;                (setf (mem-aref ,lwork-pointer :int) ,(car returned-sizes))
;;                ;; allocate and call body again
;;                (with-foreign-pointers
;;                    ,(mapcar 
;;                      (lambda (pointer atom-size returned-size)
;;                        `(,pointer (* ,atom-size ,returned-size)))
;;                      pointers atom-sizes returned-sizes)
;;                  ,@body))))))))

(defmacro with-work-area ((pointer lla-type size) &body body)
  "Allocate a work area of size lla-type elements during body,
assigning the pointer to pointer."
  (check-type pointer symbol*)
  `(with-foreign-pointer (,pointer 
			  (* ,size (foreign-size* ,lla-type)))
     ,@body))

(define-with-multiple-bindings with-work-area)

;;;; miscellaneous utility functions.

;; (defun zip-eigen% (e-val e-vec check-real? &optional
;;                    (real-type (array-lla-type e-val)))
;;   "Some LAPACK routines (yes, xGEEV, I am talking to you) return real and complex
;; parts of eigenvalues in two separate vectors.  LLA routines will try to collect these
;; in a single vector that is the concatenation of vectors containing the real and
;; imaginary parts, and this function will make a single vector out of them, given as
;; E-VAL.  If E-VEC is non-nil, it is taken to be a corresponding eigenvector matrix,
;; and eigenvectors are collected from it.  If CHECK-REAL?, eigenvectors are returned as
;; real numbers (instead of complex) when all of them are real, otherwise they are
;; always complex.  The second value returned is the matrix of eigenvectors."
;;   ;; !!! could do with more optimization and benchmarking
;;   (declare (optimize speed (safety 0)))
;;   (macrolet ((zip-eigen%% (real-type)
;;                (let ((complex-type (complex-lla-type real-type)))
;;                  `(let* ((n (aif (divides? (length e-vec) 2) it
;;                                  (error "E-VEC has to contain an even number ~
;;                                          of elements")))
;;                          (val-zipped (lla-array n ,complex-type))
;;                          (vec-zipped (when e-vec
;;                                        (let ((n^2 (expt n 2)))
;;                                          (assert (= (length e-vec) n^2))
;;                                          (lla-array n^2 ,complex-type)))))
;;                     (declare (fixnum n)
;;                              (type ,(lla-vector-type real-type) e-val e-vec)
;;                              (type ,(lla-vector-type complex-type) val-zipped))
;;                     (iter
;;                       (with index := 0)
;;                       (with all-real? := t)
;;                       (declare (iterate:declare-variables))
;;                       (let* ((realpart (aref e-val index))
;;                              (imagpart (aref e-val (+ index n)))
;;                              (real? (zerop imagpart)))
;;                         (if real?
;;                             (progn
;;                               ;; copy eigenvalue
;;                               (setf (aref val-zipped index)
;;                                     (complex realpart 0))
;;                               ;; copy a single real column
;;                               (when e-vec
;;                                 (copy-elements e-vec (cm-index2 n 0 index)
;;                                                vec-zipped (cm-index2 n 0 index) n))
;;                               ;; increment index
;;                               (incf index))
;;                             (progn
;;                               ;; copy eigenvalue and conjugate
;;                               (setf (aref val-zipped index)
;;                                     (complex realpart imagpart)
;;                                     (aref val-zipped (1+ index))
;;                                     (complex realpart (- imagpart)))
;;                               ;; copy two compex columns, conjugates of each other
;;                               (let* ((vec-index-real (cm-index2 n 0 index))
;;                                      (vec-index-imag (+ vec-index-real n)))
;; ;                                (declare (fixnum vec-index-real vec-index-imag))
;;                                 (iter
;;                                   (for (the fixnum real-index)
;;                                        :from vec-index-real
;;                                        :to vec-index-imag)
;;                                   (for (the fixnum imag-index) :from vec-index-imag)
;;                                   (declare (iterate:declare-variables))
;;                                   (let ((realpart (aref e-vec real-index))
;;                                         (imagpart (aref e-vec imag-index)))
;;                                     (setf (aref vec-zipped real-index)
;;                                           (complex realpart imagpart)
;;                                           (aref vec-zipped imag-index)
;;                                           (complex realpart (- imagpart))))))
;;                               ;; increment index
;;                               (incf index 2)))
;;                         ;; termination
;;                         (while (< index n))
;;                         (when (and all-real? (not real?))
;;                           (setf all-real? nil))
;;                         (finally
;;                          (flet ((conv (vector)
;;                                   (if (and check-real? all-real?)
;;                                       (map (lla-vector-type real-type) #'realpart vector)
;;                                       vector)))
;;                            (return (values (conv val-zipped)
;;                                            (when vec-zipped
;;                                              (make-matrix% n n 
;;                                                            (conv vec-zipped)))))))))))))
;;     (ecase real-type
;;       (:single (zip-eigen%% :single))
;;       (:double (zip-eigen%% :double)))))

;; (defun zip-eigenvalues (val-pointer n real-type complex-type
;;                         check-real-p)
;;   "Return the complex numbers stored at VAL-POINTER (N real parts,
;; followed by N imaginary parts) as a NUMERIC-VECTOR (either
;; SINGLE/DOUBLE or COMPLEX-SINGLE/COMPLEX-DOUBLE).  If CHECK-REAL-P,
;; then check if the imaginary part is 0 and if so, return a
;; NUMERIC-VECTOR-SINGLE/DOUBLE, otherwise always return a complex one.
;; The second value is non-nil if there are complex eigenvalues.  *Usage
;; note:* some LAPACK routines return real and imaginary parts of vectors
;; separately, we have to assemble them. *NOT EXPORTED*."
;;   (let ((real-p (and check-real-p 
;;                      (iter
;;                        (for i :from 0 :below n)
;;                        (always (zerop (mem-aref* val-pointer 
;;                                                  real-type 
;;                                                  (+ n i))))))))
;;     (if real-p
;;         ;; no complex eigenvalues
;;         (let ((elements (lla-array n real-type)))
;;           (iter
;;             (for i :from 0 :below n)
;;             (setf (aref elements i) 
;;                   (mem-aref* val-pointer real-type i)))
;;           (values elements nil))
;;         ;; complex eigenvalues
;;         (let ((elements (lla-array n complex-type)))
;;           (iter
;;             (for i :from 0 :below n)
;;             (setf (aref elements i) 
;;                   (complex (mem-aref* val-pointer real-type i)
;;                            (mem-aref* val-pointer real-type (+ n i)))))
;;           (values elements t)))))

;; (defun zip-eigenvectors (val-pointer vec-pointer n
;;                          real-type complex-type)
;;   "Collect complex eigenvectors from S/DGEEV.  Should only be called
;; when the second value returned by ZIP-EIGENVALUES is non-nil."
;;   (prog ((vec (lla-array (* n n) complex-type))
;;          (i 0))
;;    top
;;      (let ((column-start-index (cm-index2 n 0 i)))
;;        (if (zerop (mem-aref* val-pointer real-type (+ n i)))
;;            ;; real
;;            (progn
;;              (iter
;;                (for j :from 0 :below n)
;;                (for vec-index :from column-start-index)
;;                (setf (aref vec vec-index)
;;                      (complex 
;;                       (mem-aref* vec-pointer real-type vec-index))))
;;              (incf i))
;;            ;; complex, assemble from real +- imaginary columns
;;            (progn
;;              (iter
;;                (for j :from 0 :below n)
;;                (for vec-index :from column-start-index)
;;                (for vec-index2 :from (+ column-start-index n))
;;                (with realpart := 
;;                      (mem-aref* vec-pointer real-type vec-index))
;;                (with imagpart := 
;;                      (mem-aref* vec-pointer real-type vec-index2))
;;                (setf (aref vec vec-index) (complex realpart imagpart)
;;                      (aref vec vec-index2) (complex realpart 
;;                                                     (- imagpart))))
;;              (incf i 2))))
;;      (if (< i n)
;;          (go top)
;;          (return vec))))

;; ;;; Collecting the matrix/vector at the end.

;; (defun sum-last-rows (elements m nrhs n)
;;   "Sum & return (as a NUMERIC-VECTOR) the last (- M N) rows of an M x
;; NRHS matrix, given as a Lisp vector in column-major view.  NOTE:
;; needed to interface to LAPACK routines like xGELS."
;;   (bind ((row-sums (lla-array nrhs (real-lla-type (array-lla-type elements)))))
;;     (dotimes (col nrhs)
;;       (setf (aref row-sums col)
;;             (sse-elements% elements
;;                            (cm-index2 m n col)
;;                            (cm-index2 m m col))))
;;     row-sums))

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

(defclass lb-call-information ()
  ((procedure :documentation "see lb-call-expansion")
   (bindings :initform nil :documentation "Filtered bindings list.")
   (arrays :initform nil)
   (outputs :initform nil)
   (work-areas :initform nil)
   (preamble :initform nil)))

(defun lb-call-expansion (procedure arguments)
  "Return an expansion for calling PROCEDURE with ARGUMENTS.  Handles matrix
layout arguments (which should be preprocessed by MATRIX-LAYOUT as they are
passed directly to the function call) and error checking.  The latter works as
follows: when STATUS is given, the return value of the function is assigned to
the value named by it, and it is checked for being nonzero.  In case of an
argument error, a LAPACK-INVALID-ARGUMENT condition is raised, otherwise
CONDITION is used, with INFO."
  (let+ (((function &key layout status condition) procedure)
         (arguments (if layout
                        (cons layout arguments)
                        arguments))
         (call `(funcall ,function ,@arguments)))
    (check-type status symbol)
    (if status
        (with-unique-names (position)
          `(let ((,status ,call))
             (unless (zerop ,status)
               (if (minusp ,status)
                   (let ((,position (- ,status)))
                     (error 'lapack-invalid-argument
                            :lapack-procedure ,function
                            :argument (nth (1- ,position) ',arguments)
                            :position ,position))
                   (error ',condition :info ,status)))))
        call)))

(defun lb-process-binding (binding-form bindings)
  "Process a binding form, putting the result in bindings."
  (let+ (((variable-form &rest value-forms) binding-form))
    (if (atom variable-form)
        (push `(,variable-form ,@value-forms)
              (slot-value bindings 'bindings))
        (lb-process-binding% (car variable-form)
                             (cdr variable-form)
                             value-forms
                             bindings))))

(defmacro lb-collect ((call-information) &body slot-value-pairs)
  "Helper function for use inside lb-process-binding% methods.

   Example:

     (lb-collect (bindings)
       bindings `(,name ,value)
       arrays ...
       process `((:atom ...) value))"
  (assert (divides? (length slot-value-pairs) 2))
  (once-only (call-information)
    `(progn
       ,@(iter
           (for (slot value &rest rest) :on slot-value-pairs
                :by #'cddr)
           (collecting 
             (if (eq slot 'process)
                 `(lb-process-binding ,value ,call-information)
                 `(push ,value (slot-value ,call-information ',slot))))))))

(defun set-slot-once (call-information slot-name value)
  (assert (not (slot-boundp call-information slot-name)) ()
          "Slot ~A can only be set once." slot-name)
  (setf (slot-value call-information slot-name) value))

(defgeneric lb-process-binding% (keyword specification value-forms
                                         call-information)
  (:documentation "Process specification and add to CALL-INFORMATION."))

(defmethod lb-process-binding% (keyword specification value-forms
                                call-information)
  ;; default: process as is
  (lb-collect (call-information)
    bindings `((,keyword ,@specification) ,@value-forms)))

(defun lb-set-procedure (call-information library lla-type procedure-name 
                         layout status condition)
  (let ((procedure-name-variable (gensym* '#:procedure))
        (layout (matrix-layout library layout)))
    (set-slot-once call-information 'procedure 
                   `(,procedure-name-variable :layout ,layout :status ,status
                                              :condition ,condition))
    (lb-collect (call-information)
      bindings `(,procedure-name-variable
                 (procedure-name ,library ,lla-type
                     ,@(ensure-list procedure-name))))))

(defmethod lb-process-binding% ((keyword (eql :lapack)) specification
                                value-forms call-information)
  (let+ (((procedure-name lla-type 
              &key (layout :row-major) (status (gensym* '#:status))
              (condition 'lapack-failure))
          specification))
    (assert (null value-forms))
    (lb-set-procedure call-information :lapack lla-type procedure-name layout
                      status condition)))

(defmethod lb-process-binding% ((keyword (eql :blas)) specification
                                value-forms call-information)
  (let+ (((procedure-name lla-type &key (layout :row-major) status condition)
          specification))
    (assert (null value-forms))
    (lb-set-procedure call-information :blas lla-type procedure-name layout
                      status condition)))

;; (defmethod lb-process-binding% ((keyword (eql :atom))
;;                                 specification value-forms bindings)
;;   ;; syntax: ((:atom pointer-name lla-type) value)
;;   (bind (((pointer-name lla-type-value) specification)
;;          (lla-type-name (gensym* '#:lla-type- pointer-name))
;;          (value-name (gensym* '#:atom- pointer-name)))
;;     (lb-collect (bindings)
;;       bindings `(,value-name ,@value-forms)
;;       bindings `(,lla-type-name ,lla-type-value)
;;       atoms `(,pointer-name ,lla-type-name ,value-name))))

(defmethod lb-process-binding% ((keyword (eql :array))
                                specification value-forms call-information)
  ;; syntax: ((:matrix pointer lla-type dimensions &key
  ;;             (set-restricted? t) output) matrix)
  (let+ (((pointer-name lla-type-value &key dimensions output 
                        output-dimensions rank)
          specification)
         (lla-type-name (gensym* '#:lla-type- pointer-name))
         ((value-form) value-forms)
         (value-name (gensym* '#:value- pointer-name))
         (output-dimensions-name
          (gensym* '#:output-dimensions- pointer-name)))
    (lb-collect (call-information)
      bindings `(,value-name ,value-form)
      bindings `(,lla-type-name ,lla-type-value))
    (when output-dimensions
      (lb-collect (call-information)
        bindings `(,output-dimensions-name ,output-dimensions)))
    (lb-collect (call-information)
      arrays `(,value-name ,pointer-name ,lla-type-name :output ,output
                           ,@(when output-dimensions
                               `(:output-dimensions ,output-dimensions-name))))
    (when dimensions
      (lb-collect (call-information)
        bindings `(,(ensure-list dimensions) (array-dimensions ,value-name))))
    (when rank
      (lb-collect (call-information)
        preamble `(assert (= (array-rank ,value-name) ,rank) ()
                          "Rank of ~A is not ~A." ',value-form ,rank)))
    (awhen (and (not (eq output :copy)) output)
      (check-type it symbol*)
      (lb-collect (call-information) bindings it))))

;; (defmethod lb-process-binding% ((keyword (eql :vector))
;;                                 specification value-forms bindings)
;;   ;; syntax: ((:vector pointer-name lla-type length-spec 
;;   ;;   &optional copy-or-output) vector)
;;   (bind (((pointer-name lla-type-value length-spec &optional
;;                         output) specification)
;;          (lla-type-name (gensym* '#:lla-type- pointer-name))
;;          (value-name (gensym* '#:vector- pointer-name)))
;;     (check-type pointer-name symbol*)
;;     (check-type output (or null (eql :copy) symbol*))
;;     (lb-collect (bindings)
;;       bindings `(,value-name ,@value-forms)
;;       bindings `(,lla-type-name ,lla-type-value)
;;       arrays `(,value-name ,pointer-name ,lla-type-name :output ,output)
;;       process `((:dimension ,@(mklist length-spec))
;;                 (length ,value-name)))
;;     (awhen (and (not (eq output :copy)) output)
;;       (check-type it symbol*)
;;       (lb-collect (bindings) bindings it))))

(defmethod lb-process-binding% ((keyword (eql :output))
                                specification value-forms call-information)
  ;; syntax: ((:output pointer-name lla-type output) dimensions)
  (let+ (((pointer-name lla-type-value output-name)
          specification)
         (lla-type-name (gensym* '#:lla-type- pointer-name))
         (dimensions-name (gensym* '#:dimensions- pointer-name)))
    (check-type pointer-name symbol*)
    (check-type output-name symbol*)
    (lb-collect (call-information)
      bindings `(,dimensions-name ,@value-forms)
      bindings `(,lla-type-name ,lla-type-value)
      outputs `(,output-name ,pointer-name ,lla-type-name ,dimensions-name)
      bindings output-name)))

;; (defmethod lb-process-binding% ((keyword (eql :dimension))
;;                                 specification value-forms bindings)  
;;   ;; syntax: ((:dimension &optional var pointer) value)
;;   ;; dim-spec: NIL (not processed), VAR, (VAR), or (VAR POINTER).
;;   (when specification
;;     (bind (((variable-name &optional pointer-name) specification))
;;       (lb-collect (bindings)
;;         bindings `(,variable-name ,@value-forms))
;;       (when pointer-name
;;         (lb-collect (bindings)
;;           atoms `(,pointer-name :integer ,variable-name))))))


;; (defmethod lb-process-binding% ((keyword (eql :matrix))
;;                                 specification value-forms bindings)
;;   ;; syntax: ((:matrix pointer lla-type nrow-spec ncol-spec &key
;;   ;;             (set-restricted? t) output) matrix)
;;   ;; nrow/ncol-spec: see (:dimension ...)
;;   (bind (((pointer-name lla-type-value nrow-spec ncol-spec
;;                         &key (set-restricted? t)
;;                         output) specification)
;;          (value-name (gensym* '#:value- pointer-name)))
;;     (lb-collect (bindings)
;;       bindings `(,value-name ,@value-forms)
;;       process `((:vector ,pointer-name ,lla-type-value nil
;;                          ,output) (elements ,value-name))
;;       process `((:dimension ,@(mklist nrow-spec))
;;                 (nrow ,value-name))
;;       process `((:dimension ,@(mklist ncol-spec))
;;                 (ncol ,value-name)))
;;     (when set-restricted?
;;       (lb-collect (bindings)
;;         set-restricted? `(,value-name ,set-restricted?)))))

;; (defmethod lb-process-binding% ((keyword (eql :check))
;;                                 specification value-forms bindings)
;;   ;; syntax: ((:check info-pointer &optional (condition lapack-failure)))
;;   (bind (((:slots info) bindings)
;;          ((info-pointer &optional (condition 'lapack-failure))
;;           specification))
;;     (assert (not info) () ":CHECK already specified.")
;;     (check-type info-pointer symbol*)
;;     (check-type value-forms null)
;;     (lb-collect (bindings)
;;       info `(,info-pointer ,condition))))

;; (defmethod lb-process-binding% ((keyword (eql :work-queries))
;;                                 specification value-forms bindings)
;;   ;; syntax: ((:work-queries lwork% (work% type) ...)), where lwork% is the
;;   ;; pointer to the size of the work area, work% to the work area, and type is
;;   ;; their LLA type. !!! specifications is not expanded, just don't use any side
;;   ;; effects or depend on the order.
;;   (bind (((:slots work-queries) bindings))
;;     (assert (not work-queries) () ":WORK-QUERIES already specified.")
;;     (assert specification () ":WORK-QUERIES needs a non-empty specification.")
;;     (assert (not value-forms) () ":WORK-QUERIES does not take value forms.")
;;     (setf work-queries specification)))

(defmethod lb-process-binding% ((keyword (eql :work))
                                specification value-forms call-information)
  ;; syntax: ((:work pointer lla-type) size)
  (let+ (((pointer-name lla-type-value) specification)
         (lla-type-name (gensym* '#:lla-type- pointer-name))
         (size-name (gensym* '#:lla-type- pointer-name)))
    (check-type pointer-name symbol*)
    (lb-collect (call-information)
      bindings `(,lla-type-name ,lla-type-value)
      bindings `(,size-name ,@value-forms)
      work-areas `(,pointer-name ,lla-type-name ,size-name))))


(defmacro lb-call (binding-forms &body body)
  "The BINDING-FORMs below are captured, the rest are passed to BIND
as is.  See the comments in the corresponding LB-PROCESS-BINDING%
method.

!! todo rewrite doc here when done with reorganization

"
  (let+ ((call-information (make-instance 'lb-call-information))
         ((&flet expand (slot-name)
            (nreverse (slot-value call-information slot-name)))))
    (dolist (binding-form binding-forms)
      (lb-process-binding binding-form call-information))
    `(let+ ,(nreverse (slot-value call-information 'bindings))
       ;; define CALL macro
       (macrolet ((call (&rest arguments)
                    `(with-pinned-arrays ,',(expand 'arrays)
                       (with-array-outputs ,',(expand 'outputs)
                         (with-work-areas ,',(expand 'work-areas)
                           ,(lb-call-expansion
                             ',(slot-value call-information 'procedure)
                             arguments))))))
         ,@(expand 'preamble)
         ,@body))))
