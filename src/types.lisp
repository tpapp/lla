(in-package :lla)

;;;; Type for dimension.

(deftype dimension ()
   "Type for vector/matrix dimensions, basically a nonnegative fixnum."
  '(integer 0 #.most-positive-fixnum))

(defun check-index (index dimension)
  "Error if index is outside dimension."
  ;; ?? should this be a macro, could this signal more information
  (unless (and (<= 0 index) (< index dimension))
    (error "index ~a is outside [0,~a)" index dimension)))


;;;; Element types
;;;
;;; Some types are symbols (eg 'DOUBLE-FLOAT), some are lists (eg
;;; '(COMPLEX SINGLE-FLOAT).  CFFI uses its own naming scheme (and
;;; does not have complex types, as of Sep 2009), so we will introduce
;;; our own names -- which are sometimes appended to class and or
;;; function names --- and mappings, to and from the other naming
;;; schemes (currently only CL, as CFFI stuff is handled inside
;;; mem-aref*).

;;; NOTE Types are hardwired, because I don't think I will need more.
;;; This is not as nice/robust as LISP-MATRIX, but lookup tables
;;; would not help me much, as (1) sometimes I need to handle special
;;; cases, eg complex types with mem-aref* and (2) I would need to
;;; define a very complex DSL for possible coercions. -- Tamas

(eval-when (:compile-toplevel :load-toplevel :execute)
  (define-symbol-macro lla-types-list   ; NOT EXPORTED
      '(:single :double :complex-single :complex-double :integer))
  (deftype lla-type ()
    "All LLA types."
    `(member ,@lla-types-list))
  (defun lla-complex-p (lla-type)
    "Non-nil iff complex type."
    (or (eq lla-type :complex-single) (eq lla-type :complex-double)))
  (defun lla-double-p (lla-type)
    "Non-nil iff double-precision type."
    (or (eq lla-type :double) (eq lla-type :complex-double)))
  (defun coercible-p (lla-source-type lla-target-type)
    "Permitted coercions for LLA types.  It is guaranteed that
CL:COERCE can perform these coercions on the corresponding lisp
types. Basic summary: (1) integers can be coerced to anything, \(2)
single<->double precision coercions are possible both ways, (3) real
floats can be upgraded to complex."
    (check-type lla-source-type lla-type)
    (check-type lla-target-type lla-type)
    (cond
      ;; always valid
      ((eq lla-source-type :integer) t)
      ;; nothing else can be converted to integer
      ((eq lla-target-type :integer) nil)
      ;; no complex->real
      ((and (lla-complex-p lla-source-type)
	    (not (lla-complex-p lla-target-type)))
       nil)
      ;; all the rest should be possible
      (t t)))
  (defun coercible-pairs-list ()
      ;;   "Generate the list of all LLA (source target) pairs for which
      ;; coercible-p holds.  For internal use only, NOT EXPORTED."
      (let ((lla-types lla-types-list)
	    coercible)
	(dolist (source lla-types)
	  (dolist (target lla-types)
	    (when (coercible-p source target)
	      (push (list source target) coercible))))
	;; reverse only for cosmetic purposes
	(nreverse coercible))))

(defun lisp-type->lla-type (lisp-type &optional value-if-not-recognized)
  (cond
    ((subtypep lisp-type 'single-float) :single)
    ((subtypep lisp-type 'double-float) :double)
    ((subtypep lisp-type '(complex single-float)) :complex-single)
    ((subtypep lisp-type '(complex double-float)) :complex-double)
    ((subtypep lisp-type '(signed-byte 32)) :integer)
    (t (if (eq value-if-not-recognized 'error)
           (error "~a is not recognized as corresponding to a valid ~
  LLA type" lisp-type)
           value-if-not-recognized))))

(defun lla-type->lisp-type (lla-type)
  (case lla-type
    (:single 'single-float)
    (:double 'double-float)
    (:complex-single '(complex single-float))
    (:complex-double '(complex double-float))
    (:integer '(signed-byte 32))
    ;; !! should define & use conditions -- Tamas
    (otherwise (error "~a is not a valid LLA type" lla-type))))

(defun coerce* (value lla-type)
  "Coerce VALUE to type given by LLA-TYPE."
  (coerce value (lla-type->lisp-type lla-type)))

(defun smallest-common-target-type (lla-type-list)
  "Find the smallest supertype that all types can be coerced to.  Note
that this is always a float type because of LAPACK."
  (let ((double-p (some #'lla-double-p lla-type-list))
	(complex-p (some #'lla-complex-p lla-type-list)))
    (cond
      ((and double-p complex-p) :complex-double)
      (complex-p :complex-single)
      (double-p :double)
      (t :single))))

(defun lla-type-classifier (sequence)
  "Finds the smallest LLA-TYPE that can accomodate the elements of
  sequence.  If no such LLA-TYPE can be found, return nil."
  (ecase (reduce #'logior sequence :key
                 (lambda (x) (typecase x ;; bits: float-p, double-p, complex-p
                               ((signed-byte 32) #b000)
                               ((or single-float rational) #b100)
                               (double-float #b110)
                               ((complex single-float) #b101)
                               ((complex double-float) #b111)
                               (t (return-from lla-type-classifier nil)))))
    (#b000 :integer)
    (#b100 :single)
    (#b110 :double) 
    (#b101 :complex-single)
    (#b111 :complex-double)))     ; anything else should be impossible
