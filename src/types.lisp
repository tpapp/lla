(in-package :lla)

;;;; ** Element types
;;;;
;;;; Some types are symbols (eg 'DOUBLE-FLOAT), some are lists (eg
;;;; '(COMPLEX SINGLE-FLOAT).  CFFI uses its own naming scheme (and
;;;; does not have complex types, as of Sep 2009), so we will
;;;; introduce our own names -- which are sometimes appended to class
;;;; and or function names --- and mappings, to and from the other
;;;; naming schemes (currently only CL, as CFFI stuff is handled
;;;; inside mem-aref*).

;;;; NOTE Types are hardwired, because I don't think I will need more.
;;;; This is not as nice/robust as LISP-MATRIX, but lookup tables
;;;; would not help me much, as (1) sometimes I need to handle special
;;;; cases, eg complex types with mem-aref* and (2) I would need to
;;;; define a very complex DSL for possible coercions. -- Tamas

(eval-when (:compile-toplevel :load-toplevel :execute)
  (define-symbol-macro lla-types-list 
      ;;  Return a list of LLA types.  For internal use only, NOT EXPORTED."
      '(:single :double :complex-single :complex-double :integer)))

(eval-when (:compile-toplevel :load-toplevel :execute)
  (deftype lla-type ()
    "All LLA types."
    `(member ,@'#.lla-types-list))

  ;; ?? maybe don't need to define the types below, -p functions should
  ;; be enough -- Tamas
  
  (deftype lla-complex-type ()
    "All LLA complex (float) types."
    '(member :complex-single :complex-double))
  
  (defun lla-complex-p (lla-type)
    "Non-nil iff complex type."
    (check-type lla-type lla-type)
    (typep lla-type 'lla-complex-type))
  
  (deftype lla-double-type ()
    "All LLA double float types (real and complex)."
    '(member :double :complex-double))
  
  (defun lla-double-p (lla-type)
    "Non-nil iff double-precision type."
    (check-type lla-type lla-type)
    (typep lla-type 'lla-double-type)))


(eval-when (:compile-toplevel :load-toplevel :execute)
  (defun coercible-p (lla-source-type lla-target-type)
    "Permitted coercions for LLA types.  It is guaranteed that CL:COERCE
can perform these coercions on the corresponding lisp types.

Basic summary:

- integers can be coerced to anything
- single<->double precision coercions are possible both ways
- real floats can be upgraded to complex."
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


(defmacro for-coercible-pairs ((source target) form)
  "Replace source and target with corresponding coercible pairs
in form for each coercible pair in turn, wrap the result in a progn."
  (cons 'progn
	(mapcar (lambda (pair)
		  (sublis (list (cons source (first pair))
				(cons target (second pair)))
			  form))
		(coercible-pairs-list))))

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
  "Finds the smallest lla type that can accomodate the elements of
  list (after being converted with coerce).  If no such lla-type can
  be found, return nil."
  (let ((typecode (reduce (lambda (x y)
                            ;; parallel maximum on a cons
                            (cons (max (car x) (car y)) (max (cdr x) (cdr y))))
                          sequence
                          :key (lambda (x)
                                 ;; conses: precision (0: integer, 1:
                                 ;; single, 2: double) and complex indicator
                                 (cond
                                   ((typep x '(signed-byte 32)) '(0 . 0))
                                   ((typep x '(or single-float rational)) '(1 . 0))
                                   ((typep x 'double-float) '(2 . 0))
                                   ((typep x '(complex single-float)) '(1 . 1))
                                   ((typep x '(complex double-float)) '(2 . 1))
                                   (t (return-from lla-type-classifier nil)))))))
    (cond
      ((equal typecode '(0 . 0)) :integer)
      ((equal typecode '(1 . 0)) :single)
      ((equal typecode '(2 . 0)) :double) 
      ((equal typecode '(1 . 1)) :complex-single)
      ((equal typecode '(2 . 1)) :complex-double)
      (t (error "internal error, this should never happen")))))

(defun to-lla-type (number lla-type)
  "Return number coerced to lla-type."
  (coerce number (lla-type->lisp-type lla-type)))
