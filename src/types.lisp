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
  (defparameter +lla-type-list+
    '(:single :double :complex-single :complex-double :integer)
    "This should be a constant, if CL supported lists as constants.
This is not exported: types are hardwired into LLA, they are just
enumerated here for convenience."))

(deftype lla-type ()
  "All LLA types."
  `(member ,@+lla-type-list+))

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
  (let ((lla-types +lla-type-list+)
        coercible)
    (dolist (source lla-types)
      (dolist (target lla-types)
        (when (coercible-p source target)
          (push (list source target) coercible))))
    ;; reverse only for cosmetic purposes
    (nreverse coercible)))

(defgeneric lla-type (numeric-vector)
  (:documentation "Return the lla-type of the elements of the
  object."))

(define-condition not-within-lla-type (error)
  ()
  (:documentation "Could not classify given type as a subtype of an LLA-TYPE."))

(defun lisp-type->lla-type (lisp-type &optional value-if-not-recognized)
  (cond
    ((subtypep lisp-type 'single-float) :single)
    ((subtypep lisp-type 'double-float) :double)
    ((subtypep lisp-type '(complex single-float)) :complex-single)
    ((subtypep lisp-type '(complex double-float)) :complex-double)
    ((subtypep lisp-type '(signed-byte 32)) :integer)
    (t (if (eq value-if-not-recognized 'error)
           (error 'not-within-lla-type)
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

(defmacro append-lla-type (prefix lla-type)
  "Return prefix-LLA-TYPE."
  `(ecase ,lla-type
     (:integer ',(make-symbol* prefix '-integer))
     (:single ',(make-symbol* prefix '-single))
     (:double ',(make-symbol* prefix '-double))
     (:complex-single ',(make-symbol* prefix '-complex-single))
     (:complex-double ',(make-symbol* prefix '-complex-double))))

;;;; automatic type classification
;;;
;;; We find the smallest common target type, the type we can all
;;; coerce to, using binary or.  Bits: float-p, double-p, complex-p

(defconstant +integer+ #b000)
(defconstant +single+ #b100)
(defconstant +double+ #b110)
(defconstant +complex-single+ #b101)
(defconstant +complex-double+ #b111)

(defconstant +forced-float+ :double)

(declaim (inline lla-type->binary-code binary-code->lla-type
                 common-target-type))

(defun lla-type->binary-code (type)
  "Convert LLA type to bits."
  (ecase type
    (:integer +integer+)
    (:single +single+)
    (:double +double+)
    (:complex-single +complex-single+)
    (:complex-double +complex-double+)))

(defvar *force-float* t "If non-nil, sequences with automatically
detected :integer types will be converted to +FORCED-FLOAT+.  This
should only influence type autodetection, not operations.")

(defun binary-code->lla-type (binary-code &optional force-float-p)
  "Convert bits to LLA type."
  (ecase binary-code
    (#.+integer+ (if force-float-p +forced-float+ :integer))
    (#.+single+ :single)
    (#.+double+ :double)
    (#.+complex-single+ :complex-single)
    (#.+complex-double+ :complex-double)))

(defun common-target-type (&rest types)
  "Find the smallest supertype that all types can be coerced to.  Does
NOT force floats."
  (binary-code->lla-type
   (reduce #'logior types :key #'lla-type->binary-code)))

(defun find-element-type (sequence &optional (force-float-p *force-float*))
  "Finds the smallest LLA-TYPE that can accomodate the elements of
  sequence.  If no such LLA-TYPE can be found, return nil."
  (binary-code->lla-type
   (reduce #'logior sequence :key (compose #'lla-type->binary-code #'lisp-type->lla-type #'type-of))
   force-float-p))

(defun infer-lla-type (lla-type initial-contents &optional (force-float-p *force-float*))
  "Infer LLA-TYPE from given type and/or INITIAL-CONTENTS.  Useful for
MAKE-NV and MAKE-MATRIX."
  (cond
    (lla-type lla-type)
    ((null initial-contents) (error "Could not infer LLA-TYPE without initial-contents."))
    ((numberp initial-contents) (let ((lla-type (lisp-type->lla-type (type-of initial-contents))))
                                  (if (and force-float-p (eq lla-type :integer))
                                      +forced-float+ lla-type)))
    ((typep initial-contents 'sequence) (find-element-type initial-contents force-float-p))
    (t (error "~A is not valid as INITIAL-CONTENTS." initial-contents))))


;;;; macros for defining subclasses of NUMERIC-VECTOR-*


(defmacro expand-for-lla-types ((typevar &key (prologue '(progn))
                                         (exclude-integer-p nil))
                                &body form)
  "Expand FORM (using EVAL) with TYPEVAR bound to all possible LLA
types, return the results inside a (,@PROLOGUE ...), PROGN by
default."
  (when (cadr form)
    ;; FORM is a &body argument only for saner indentation.
    (error "Multiple forms provided."))
  `(,@prologue
    ,@(mapcar (lambda (typename)
                (eval `(let ((,typevar ',typename)) ,(car form))))
              (if exclude-integer-p
                  '(:single :double :complex-single :complex-double)
                  +lla-type-list+))))

(defmacro define-lla-class (class &optional (superclasses (list class)))
  "Define subclasses of all NUMERIC-VECTOR types and given
superclasses, with appropriate names.  If no superclasses are given,
class is used instead."
  `(progn ,@(mapcar 
             (lambda (lla-type)
               `(defclass ,(make-symbol* class "-" lla-type)
                      (,@superclasses ,(append-lla-type numeric-vector lla-type))
                  ()))
             +lla-type-list+)))
