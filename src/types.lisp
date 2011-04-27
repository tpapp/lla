(in-package :lla)

;; (define-condition not-within-lla-type (error)
;;   ()
;;   (:documentation "Could not classify given type as a subtype of an
;;   LLA-TYPE."))

(define-condition invalid-lla-type (error)
  ()
  (:documentation "The given type is not a valid LLA type."))

(eval-when (:compile-toplevel :load-toplevel :execute)
  (defun lla-types ()
    "Return a list of valid LLA types."
    '(:integer :single :double :complex-single :complex-double))
  (defun lla-to-lisp-type (lla-type)
    "Return the lisp type spec corresponding to LLA-TYPE.  Signals an error for
invalid LLA types."
    (case lla-type
      (:single 'single-float)
      (:double 'double-float)
      (:complex-single '(complex single-float))
      (:complex-double '(complex double-float))
      (:integer '(signed-byte #-lla-int64 32 #+lla-int64 64))
      (otherwise (error 'invalid-lla-type))))
  (defun atom-lla-type (atom)
    "If atom is a subtype of a Lisp type corresponding to an LLA type, return that,
otherwise NIL."
    (typecase atom
      (single-float :single)
      (double-float :double)
      ((complex single-float) :complex-single)
      ((complex double-float) :complex-double)
      ((signed-byte #-lla-int64 32 #+lla-int64 64) :integer))))

(defun lla-to-lisp-type* (lla-or-lisp-type)
  "Like LLA-TO-LISP-TYPE, but pass through invalid LLA types unchanged."
  (handler-case (lla-to-lisp-type lla-or-lisp-type)
    (invalid-lla-type () lla-or-lisp-type)))

(defun lla-array-element-type* (lla-or-lisp-type)
  "Like LLA-ARRAY-ELEMENT-TYPE, but also allow CL types."
  (handler-case (lla-array-element-type lla-or-lisp-type)
    (invalid-lla-type () (upgraded-array-element-type lla-or-lisp-type))))

;; (deftype lla-type ()
;;   "All LLA types."
;;   `(member :single :double :complex-single :complex-double :integer))

(defun lla-complex? (lla-type)
  "Non-nil iff complex type."
  (or (eq lla-type :complex-single) (eq lla-type :complex-double)))

(defun lla-double? (lla-type)
  "Non-nil iff double-precision type."
  (or (eq lla-type :double) (eq lla-type :complex-double)))

(defun real-lla-type (lla-type)
    "Return the lla-type of (* x (conjugate x), when x is of LLA-TYPE."
    (ecase lla-type
      (:integer :integer)
      (:single :single)
      (:double :double)
      (:complex-single :single)
      (:complex-double :double)
      ;; ((nil) nil)
))

(defun complex-lla-type (lla-type)
  "Return the complex type corresponding to LLA-TYPE."
  (ecase lla-type
    ((:integer :single :complex-single) :complex-single)
    ((:double :complex-double) :complex-double)))

(declaim (inline zero* coerce* epsilon*))
(defun zero* (type)
  "Return 0, coerced to the desired TYPE (extended with LLA types)."
  (case type
    (:single 0.0)
    (:double 0d0)
    (:complex-single #C(0.0 0.0))
    (:complex-double #C(0.0d0 0.0d0))
    ((:integer nil) 0)
    (otherwise (coerce 0 type))))

(defun coerce* (value lla-type)
  "Coerce VALUE to type given by LLA-TYPE, NIL implies no coercion."
  (coerce value (lla-to-lisp-type* lla-type)))

(defun epsilon* (lla-type)
  "Return the float epsilon for the given LLA-TYPE, signal an error
for :INTEGER."
  (ecase lla-type
    ((:single :complex-single) single-float-epsilon)
    ((:double :complex-double) double-float-epsilon)))

;; representable LLA types: discover properties of the implementation

(eval-when (:compile-toplevel :load-toplevel :execute)
  (defmacro define-array-functions% ()
    (let ((representable-table
           (loop
             for lla-type :in (lla-types)
             collect (let* ((lisp-type (lla-to-lisp-type lla-type))
                            (upgraded-type (upgraded-array-element-type lisp-type)))
                       (cons (and (subtypep lisp-type upgraded-type)
                                  (subtypep upgraded-type lisp-type))
                             lla-type)))))
      `(progn
         (defun lla-array-element-type (lla-type)
           "Return the array element type corresponding to LLA-TYPE.  If the latter
is not representable in the implementation, return T.  Signal an error for invalid
LLA-TYPE."
           (case lla-type
             ,@(loop
                 for (representable? . lla-type) :in representable-table
                 collect `(,lla-type ',(if representable? 
                                           (lla-to-lisp-type lla-type)
                                           t)))
             (otherwise (error 'invalid-lla-type))))
         (defun array-manifest-lla-type (array)
           "If the element-type of ARRAY corresponds to a LLA-TYPE, return that,
otherwise NIL."
           (let ((element-type (array-element-type array)))
             (cond
               ,@(loop
                   for (representable? . lla-type) :in representable-table
                   when representable?
                     collect `((subtypep element-type ',(lla-to-lisp-type lla-type))
                                         ,lla-type))))))))
  (define-array-functions%))

(defun make-array* (dimensions element-type &optional initial-element)
  "Create an array with given LLA or lisp type and dimensions.  INITIAL-ELEMENT is
used to fill the array when given."
  (let ((lisp-type (lla-array-element-type* element-type)))
    (if initial-element
        (make-array dimensions :element-type lisp-type
                    :initial-element (coerce* initial-element element-type))
        (make-array dimensions :element-type lisp-type))))

(defun convert-array* (array element-type)
  "Convert an ARRAY to a given LLA type."
  (aprog1 (make-array* (array-dimensions array) element-type)
    (map-into (flatten-array it) (lambda (x) (coerce* x element-type))
              (flatten-array array))))

(defun maybe-convert-array* (array element-type copy?)
  "Convert array only when types mismatch or COPY? is non-NIL."
  (if (and (equal (array-element-type array) (lla-array-element-type* element-type))
           (not copy?))
      array
      (convert-array* array element-type)))

(defun common-lla-type (objects &key force-float? double?)
  "Find narrowest common LLA type for OBJECTS.  FORCE-FLOAT? turns integers into
floats, FORCE-DOUBLE? determines whether to return :DOUBLE or :SINGLE for rationals,
mutatis mutandis for complex rationals.  When OBJECTS do not have a common type
representabel in the LLA framework, return NIL."
  ;; implementation note: + takes care of type determination via contagion.  Crude,
  ;; but works fine and is fast.
  (let ((float-type (if double? :double :single))
        (sum (reduce #'+
                     (mapcar (lambda (object)
                               (etypecase object
                                 (array (if (array-manifest-lla-type object)
                                            (row-major-aref object 0)
                                            (reduce #'+ (flatten-array object))))
                                 (number object)))
                             objects))))
    (typecase sum
      (single-float :single)
      (double-float :double)
      (integer (if force-float? float-type :integer))
      (rational float-type)
      ((complex single-float) :complex-single)
      ((complex double-float) :complex-double)
      ((complex rational) (if double? :complex-double :complex-single)))))

(defgeneric pack (object &key force-float? double?)
  (:documentation "The purpose of this is to return objects for immediate
  detection of element types, affecting speed favorably.  Vectors are converted
  to simple-array.  The returned object may share structure with the
  original."))

(defmethod pack ((array array) &key force-float? double?)
  (let ((lla-type (common-lla-type (list array) 
                                   :force-float? force-float? :double? double?)))
    (if (or (eq (array-manifest-lla-type array) lla-type)
            (eq (lla-array-element-type lla-type) t))
        array
        (convert-array* array lla-type))))

(define-modify-macro packf () pack "PACK place.")

(declaim (inline conjugate-square))

(defun conjugate-square (number)
  "Number multiplied by its complex conjugate."
  (* (conjugate number) number))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; moving wall, clean up below this
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;; Notes on LLA's optimization model
;;;
;;; LLA is designed to allow the use of Lisp arrays (most notably vectors) almost
;;; everywhere.  All operations work directly on vectors, or lightweight wrappers
;;; around vectors.
;;;
;;; 1. Upgraded element types
;;;
;;; When interfacing with BLAS/LAPACK, LLA tries to make use of the
;;; implementation's facilities for avoiding copying into foreign
;;; memory (array "pinning"), whenever possible.  However, for
;;; deciding which BLAS/LAPACK function to call, the type of elements
;;; in vectors or matrices has to be established.
;;;
;;; For this purpose, LLA recognizes five numeric types as special, and denotes them
;;; by the following keywords: :INTEGER (the integer type used by your BLAS/LAPACK
;;; library, most likely 32-bit), :SINGLE, :DOUBLE, :COMLEX-SINGLE and
;;; :COMPLEX-DOUBLE, for single- and double-floats, either real or complex.  NIL is
;;; used for arrays which were not recognized as either of these.
;;;
;;; Some implementations upgrade these element types to themselves in
;;; arrays, in which case LLA can select types very quickly.  In case
;;; the implementation upgrades to something else (eg T, but not
;;; necessarily), LLA has to sweep the array before BLAS/LAPACK calls
;;; before to establish the common element type.
;;;
;;; 2. Simple arrays
;;;
;;; Elements in LLA objects are represented by rank-one arrays.  LLA
;;; prefers if these arrays are simple, this makes certain
;;; optimizations possible by declaring this type.  Because of this,
;;; LLA objects that encapsulate arrays always enforce the restriction
;;; that they are simple (specifically, of type simple-array1, which
;;; is a shorthand for (simple-array * (*)).  For vectors which are
;;; not encapsulated in LLA objects, either a conversion will take
;;; place silently, or an error message will result (the former is
;;; more common).
;;;
;;; 3. Optimal way of using LLA
;;;
;;; LLA will gobble any kind of vector if the elements make sense for
;;; the given operation.  But you can avoid unnecessary
;;; conversions/checks by adhering to the following guidelines:
;;;
;;;   a. Keep you vectors rank 1, SIMPLE-ARRAY, and preferably of the
;;;   narrowest element type supported by your implementation.  See
;;;   LLA-ARRAY and AS-SIMPLE-ARRAY1.
;;;
;;;   b. When using matrices, use DENSE-MATRIX-LIKE objects.  LLA is
;;;   column-major, and 
;;;
;;; Also see *lla-warn-suboptimal*.

;;; Non-nil symbols.



;; ;;; Type for dimension.

;; (deftype dimension ()
;;    "Type for vector/matrix dimensions, basically a nonnegative fixnum."
;;   '(integer 0 #.array-total-size-limit))

;; (defun check-index (index dimension)
;;   "Error if index is outside dimension."
;;   ;; ?? should this be a macro, could this signal more information
;;   (unless (within? 0 index dimension)
;;     (error "index ~a is outside [0,~a)" index dimension)))

;; ;;; warnings for suboptimal operations

;; (defvar *lla-warn-suboptimal* nil
;;   "When set to non-NIL, emit a warning when suboptimal operations are peformed.
;; If set to BREAK, break will be used to halt execution.  If NIL, no warning is
;; generated.

;; Note: these operations can usually be improved by using simple-arrays or rank 1
;; and supported LLA element types.  Note: this is intended for debugging your
;; application on implementations that support the appropriate specialized arrays.
;; On implementations that don't (eg CLISP), setting this will just emit a lot of
;; warning messages regardless of what you do.")

;; (define-condition lla-suboptimal-warning (warning)
;;   ((message :initarg :message)))

;; (defmethod print-object ((warning lla-suboptimal-warning) stream)
;;   (format stream "Suboptional LLA operation: ~A."
;;           (slot-value warning 'message)))

;; (defun lla-warn-suboptimal (message)
;;   "Signal suboptimal LLA operations, according to *LLA-WARN-SUBOPTIMAL*."
;;   (case *lla-warn-suboptimal*
;;     ((nil))
;;     (break (break "Suboptional LLA operation: ~A." message))
;;     (otherwise (warn 'lla-suboptimal-warning :message message)))
;;   (values))

;; ;;; Element types
;; ;;;
;; ;;; In CL, some types are symbols (eg 'DOUBLE-FLOAT), some are lists
;; ;;; (eg '(COMPLEX SINGLE-FLOAT).  CFFI uses its own naming scheme (and
;; ;;; does not have complex types, as of Sep 2009), so we will introduce
;; ;;; our own names -- which are sometimes appended to class and or
;; ;;; function names --- and mappings, to and from the other naming
;; ;;; schemes (currently only CL, as CFFI stuff is handled inside
;; ;;; mem-aref*).
;; ;;;
;; ;;; NOTE Types are hardwired, because I don't think I will need more.
;; ;;; This is not as nice/robust as LISP-MATRIX, but lookup tables
;; ;;; would not help me much, as (1) sometimes I need to handle special
;; ;;; cases, eg complex types with mem-aref* and (2) I would need to
;; ;;; define a very complex DSL for possible coercions. -- Tamas



;; ;; (defun atom-representable-lla-type (atom)
;; ;;   "If the (upgraded) type of atom is supported by the implementation,
;; ;; return the corresponding LLA type, otherwise NIL."
;; ;;   (typecase atom
;; ;;     #+lla-vector-integer ((signed-byte #-lla-int64 32 #+lla-int64 64) :integer)
;; ;;     #+lla-vector-single (single-float :single)
;; ;;     #+lla-vector-double (double-float :double)
;; ;;     #+lla-vector-complex-single ((complex single-float) :complex-single)
;; ;;     #+lla-vector-complex-double ((complex double-float) :complex-double)))


;; ;; (defun similar-simple-array? (first &rest rest)
;; ;;   (let ((element-type (array-element-type first)))
;; ;;     (every (lambda (array)
;; ;;              (and (simple-array? array)
;; ;;                   (equal (array-element-type array) element-type)))
;; ;;            rest)))

;; ;; (defmacro with-vector-type-expansion ((vector &key
;; ;;                                               other-vectors
;; ;;                                               (simple-test t)
;; ;;                                               (vector? t)
;; ;;                                               (lla-type (gensym "LLA-TYPE")))
;; ;;                                       body-generator)
;; ;;   "Expand based on the LLA type of vector, declaring the vector (and
;; ;; other vectors, if given) to be of this type.  Only expand w/ type
;; ;; declarations for upgraded element type is supported by the
;; ;; implementation.  Also, these cases apply only when VECTOR is a
;; ;; SIMPLE-VECTOR, and SIMPLE-TEST is satisfied.

;; ;; The body for each case is generated by converting body-generator to a
;; ;; function, and calling it on the appropriate LLA type, including NIL as
;; ;; the fallback case.

;; ;; Usage note: it is implicitly assumed that OTHER-VECTORS are
;; ;; SIMPLE-VECTORs, created with the same element type as VECTOR.  If this
;; ;; is cannot be known for certain, use SIMPLE-TEST: :OTHER-VECTORS
;; ;; implements this behavior by default."
;; ;;   (check-type vector symbol)
;; ;;   (assert (every #'symbolp other-vectors))
;; ;;   (setf body-generator (coerce body-generator 'function))
;; ;;   (with-unique-names (simple?)
;; ;;     (flet ((clause (case-lla-type)
;; ;;              `((and ,simple? (eq ,lla-type ,case-lla-type))
;; ;;                (locally 
;; ;;                    (declare (type (simple-array 
;; ;;                                    ,(lla->lisp-type case-lla-type)
;; ;;                                    ,(if vector? '(*) '*))
;; ;;                                   ,vector ,@other-vectors))
;; ;;                  ,(funcall body-generator case-lla-type)))))
;; ;;       `(let* ((,lla-type (array-lla-type ,vector))
;; ;;               (,simple? (and (simple-array? ,vector)
;; ;;                              ,(if (eq simple-test :other-vectors)
;; ;;                                   `(similar-simple-array? ,vector
;; ;;                                                           ,@other-vectors)
;; ;;                                   simple-test))))
;; ;;          (cond
;; ;;            #+lla-vector-integer
;; ;;            ,(clause :integer)
;; ;;            #+lla-vector-single
;; ;;            ,(clause :single)
;; ;;            #+lla-vector-double
;; ;;            ,(clause :double)
;; ;;            #+lla-vector-complex-single
;; ;;            ,(clause :complex-single)
;; ;;            #+lla-vector-complex-double
;; ;;            ,(clause :complex-double)
;; ;;            (t (muffle-optimization-notes
;; ;;                 (lla-warn-suboptimal "non-specialized vector")
;; ;;                 ,(funcall body-generator nil))))))))


;; ;; (defmacro with-vector-type-declarations ((vector &key
;; ;;                                                  other-vectors
;; ;;                                                  (simple-test t)
;; ;;                                                  (vector? t)
;; ;;                                                  (lla-type (gensym "LLA-TYPE")))
;; ;;                                          &body body)
;; ;;   "Like WITH-VECTOR-TYPE-EXPANSION, but constant body-generator taken
;; ;; from BODY."
;; ;;   `(with-vector-type-expansion (,vector
;; ;;                                 :other-vectors ,other-vectors
;; ;;                                 :simple-test ,simple-test
;; ;;                                 :vector? ,vector?
;; ;;                                 :lla-type ,lla-type)
;; ;;      (lambda (lla-type)
;; ;;        (declare (ignore lla-type))
;; ;;        `(progn ,',@body))))

;; ;; (defun common-lla-type (lla-type1 lla-type2)
;; ;;   "Return a common LLA type that both of the given types can be
;; ;; coerced to, or NIL if that is not possible.  Arguments have to be
;; ;; valid LLA types, which is not checked."
;; ;;   (cond
;; ;;     ((or (null lla-type1) (null lla-type2))
;; ;;      nil)
;; ;;     ((eq lla-type1 lla-type2)
;; ;;      lla-type1)
;; ;;     (t (let ((complex? (or (lla-complex? lla-type1)
;; ;;                            (lla-complex? lla-type2)))
;; ;;              (double? (or (lla-double? lla-type1)
;; ;;                           (lla-double? lla-type2))))
;; ;;          (if complex?
;; ;;              (if double? :complex-double :complex-single)
;; ;;              (if double? :double :single))))))

;; ;; (defgeneric elements (object)
;; ;;   (:documentation "Return the elements of an LLA object as an array or a scalar."))

;; ;; (defparameter *lla-force-float* nil
;; ;;   "Determines whether LLA functions convert integers to float.")

