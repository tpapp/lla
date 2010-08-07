(in-package :lla)

;;; Notes on LLA's optimization model
;;;
;;; LLA is designed to allow the use of Lisp arrays (most notably
;;; vectors) almost everywhere.  Many operations work directly on
;;; vectors, and most LLA objects are nothing but a lightweight
;;; wrapper for such vectors.
;;;
;;; 1. Upgraded element types
;;;
;;; When interfacing with BLAS/LAPACK, LLA tries to make use of the
;;; implementation's facilities for avoiding copying into foreign
;;; memory (array "pinning"), whenever possible.  However, for
;;; deciding which BLAS/LAPACK function to call, the type of elements
;;; in vectors or matrices has to be established.
;;;
;;; For this purpose, LLA recognizes five numeric types as special,
;;; and denotes them by the following keywords: :INTEGER (for 32 or
;;; 64-bit integers, depending on the platform), :SINGLE, :DOUBLE,
;;; :COMLEX-SINGLE and :COMPLEX-DOUBLE, for single- and double-floats,
;;; either real or complex.  NIL is used for arrays which were not
;;; recognized as either of these.
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
;;;   LLA-VECTOR and AS-SIMPLE-ARRAY1.
;;;
;;;   b. When using matrices, use DENSE-MATRIX-LIKE objects.  LLA is
;;;   column-major, and 
;;;
;;; Also see *lla-warn-suboptimal*.

;;; Non-nil symbols.

(deftype symbol* ()
  '(and symbol (not null)))

(defun symbolp* (object)
  (typep object 'symbol*))

;;; Type for dimension.

(deftype dimension ()
   "Type for vector/matrix dimensions, basically a nonnegative fixnum."
  '(integer 0 #.array-total-size-limit))

(defun check-index (index dimension)
  "Error if index is outside dimension."
  ;; ?? should this be a macro, could this signal more information
  (unless (within? 0 index dimension)
    (error "index ~a is outside [0,~a)" index dimension)))

;;; warnings for suboptimal operations

(defvar *lla-warn-suboptimal* nil
  "When set to non-NIL, emit a warning when suboptimal operations are peformed.
If set to BREAK, break will be used to halt execution.  If NIL, no warning is
generated.

Note: these operations can usually be improved by using simple-arrays or rank 1
and supported LLA element types.  Note: this is intended for debugging your
application on implementations that support the appropriate specialized arrays.
On implementations that don't (eg CLISP), setting this will just emit a lot of
warning messages regardless of what you do.")

(define-condition lla-suboptimal-warning (warning)
  ((message :initarg :message)))

(defmethod print-object ((warning lla-suboptimal-warning) stream)
  (format stream "Suboptional LLA operation: ~A."
          (slot-value warning 'message)))

(defun lla-warn-suboptimal (message)
  "Signal suboptimal LLA operations, according to *LLA-WARN-SUBOPTIMAL*."
  (case *lla-warn-suboptimal*
    ((nil))
    (break (break "Suboptional LLA operation: ~A." message))
    (otherwise (warn 'lla-suboptimal-warning :message message)))
  (values))

;;; Element types
;;;
;;; In CL, some types are symbols (eg 'DOUBLE-FLOAT), some are lists
;;; (eg '(COMPLEX SINGLE-FLOAT).  CFFI uses its own naming scheme (and
;;; does not have complex types, as of Sep 2009), so we will introduce
;;; our own names -- which are sometimes appended to class and or
;;; function names --- and mappings, to and from the other naming
;;; schemes (currently only CL, as CFFI stuff is handled inside
;;; mem-aref*).
;;;
;;; NOTE Types are hardwired, because I don't think I will need more.
;;; This is not as nice/robust as LISP-MATRIX, but lookup tables
;;; would not help me much, as (1) sometimes I need to handle special
;;; cases, eg complex types with mem-aref* and (2) I would need to
;;; define a very complex DSL for possible coercions. -- Tamas

(deftype lla-type ()
  "All LLA types."
  `(member :single :double :complex-single :complex-double :integer))

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
      ((nil) nil)))

(defun complex-lla-type (lla-type)
  "Return the complex type corresponding to LLA-TYPE."
  (ecase lla-type
    ((:integer :single :complex-single) :complex-single)
    ((:double :complex-double) :complex-double)
    ((nil) nil)))

(define-condition not-within-lla-type (error)
  ()
  (:documentation "Could not classify given type as a subtype of an
  LLA-TYPE."))

(define-condition invalid-lla-type (error)
  ()
  (:documentation "The given type is not a valid LLA type."))

(eval-when (:compile-toplevel :load-toplevel :execute)
  (defun lla->lisp-type (lla-type)
    "Return the lisp type spec corresponding to LLA-TYPE.  Signals an
error for invalid LLA types."
    (case lla-type
      (:single 'single-float)
      (:double 'double-float)
      (:complex-single '(complex single-float))
      (:complex-double '(complex double-float))
      (:integer '(signed-byte #-lla-int64 32 #+lla-int64 64))
      ((nil) t)
      (otherwise (error 'invalid-lla-type)))))

(defun lla-vector-type (lla-type)
  "Return the type of a vector with elements conforming to LLA-TYPE.  May be
upgraded by the implementation."
  `(simple-array ,(lla->lisp-type lla-type) (*)))

(eval-when (:compile-toplevel :load-toplevel :execute)
  (labels ((type= (type1 type2)
             "Test if two types are equal (ie subtypes of each other)."
             (and (subtypep type1 type2)
                  (subtypep type2 type1)))
           (check-upgraded (lla-type feature)
             (let ((lisp-type (lla->lisp-type lla-type)))
               (when (type= (upgraded-array-element-type lisp-type) lisp-type)
                 (pushnew feature *features*)))))
    (check-upgraded :integer :lla-vector-integer)
    (check-upgraded :single :lla-vector-single)
    (check-upgraded :double :lla-vector-double)
    (check-upgraded :complex-single :lla-vector-complex-single)
    (check-upgraded :complex-double :lla-vector-complex-double)))

(defun representable-lla-type (type)
  "If TYPE corresponds to an LLA type representable in an array, return that,
otherwise NIL."
  (cond
    #+lla-vector-integer 
    ((subtypep type '(signed-byte #-lla-int64 32 #+lla-int64 64)) :integer)
    #+lla-vector-single 
    ((subtypep type 'single-float) :single)
    #+lla-vector-double 
    ((subtypep type 'double-float) :double)
    #+lla-vector-complex-single 
    ((subtypep type '(complex single-float)) :complex-single)
    #+lla-vector-complex-double 
    ((subtypep type '(complex double-float)) :complex-double)))

(defun array-lla-type (array)
  "If the array element type corresponds to an LLA type, return that, otherwise
NIL."
  (representable-lla-type (array-element-type array)))

(defun lla-vector (length lla-type &optional initial-element)
  "Create a vector with given LLA type and length.  Initial-element will be
coerced to the appropriate type and used to fill the vector (the default is 0).
It can also be a function, in which case it will be called to fill each element,
traversing from index 0 to length-1."
  (bind ((lisp-type (lla->lisp-type lla-type))
         ((:flet make (&optional initial-element))
          (if initial-element
              (make-array length :element-type lisp-type
                          :initial-element initial-element)
              (make-array length :element-type lisp-type))))
    (typecase initial-element
      (null (make))
      (function (aprog1 (make)
                  (dotimes (index length)
                    (setf (aref it index) (funcall initial-element)))))
      (otherwise (make (coerce* initial-element lla-type))))))

(defun atom-representable-lla-type (atom)
  "If the (upgraded) type of atom is supported by the implementation,
return the corresponding LLA type, otherwise NIL."
  (typecase atom
    #+lla-vector-integer ((signed-byte #-lla-int64 32 #+lla-int64 64) :integer)
    #+lla-vector-single (single-float :single)
    #+lla-vector-double (double-float :double)
    #+lla-vector-complex-single ((complex single-float) :complex-single)
    #+lla-vector-complex-double ((complex double-float) :complex-double)))

(declaim (inline zero* coerce* epsilon*))
(defun zero* (lla-type)
  "Return 0, coerced to the desired LLA-TYPE."
  (ecase lla-type
    (:single 0.0)
    (:double 0d0)
    (:complex-single #C(0.0 0.0))
    (:complex-double #C(0.0d0 0.0d0))
    ((:integer nil) 0)))

(defun coerce* (value lla-type)
  "Coerce VALUE to type given by LLA-TYPE, NIL implies no coercion."
  (coerce value (lla->lisp-type lla-type)))

(defun epsilon* (lla-type)
  "Return the float epsilon for the given LLA-TYPE, signal an error
for :INTEGER."
  (ecase lla-type
    ((:single :complex-single) single-float-epsilon)
    ((:double :complex-double) double-float-epsilon)))

(defun similar-simple-array? (first &rest rest)
  (let ((element-type (array-element-type first)))
    (every (lambda (array)
             (and (simple-array? array)
                  (equal (array-element-type array) element-type)))
           rest)))

(defmacro with-vector-type-expansion ((vector &key
                                              other-vectors
                                              (simple-test t)
                                              (vector? t)
                                              (lla-type (gensym "LLA-TYPE")))
                                      body-generator)
  "Expand based on the LLA type of vector, declaring the vector (and
other vectors, if given) to be of this type.  Only expand w/ type
declarations for upgraded element type is supported by the
implementation.  Also, these cases apply only when VECTOR is a
SIMPLE-VECTOR, and SIMPLE-TEST is satisfied.

The body for each case is generated by converting body-generator to a
function, and calling it on the appropriate LLA type, including NIL as
the fallback case.

Usage note: it is implicitly assumed that OTHER-VECTORS are
SIMPLE-VECTORs, created with the same element type as VECTOR.  If this
is cannot be known for certain, use SIMPLE-TEST: :OTHER-VECTORS
implements this behavior by default."
  (check-type vector symbol)
  (assert (every #'symbolp other-vectors))
  (setf body-generator (coerce body-generator 'function))
  (with-unique-names (simple?)
    (flet ((clause (case-lla-type)
             `((and ,simple? (eq ,lla-type ,case-lla-type))
               (locally 
                   (declare (type (simple-array 
                                   ,(lla->lisp-type case-lla-type)
                                   ,(if vector? '(*) '*))
                                  ,vector ,@other-vectors))
                 ,(funcall body-generator case-lla-type)))))
      `(let* ((,lla-type (array-lla-type ,vector))
              (,simple? (and (simple-array? ,vector)
                             ,(if (eq simple-test :other-vectors)
                                  `(similar-simple-array? ,vector
                                                          ,@other-vectors)
                                  simple-test))))
         (cond
           #+lla-vector-integer
           ,(clause :integer)
           #+lla-vector-single
           ,(clause :single)
           #+lla-vector-double
           ,(clause :double)
           #+lla-vector-complex-single
           ,(clause :complex-single)
           #+lla-vector-complex-double
           ,(clause :complex-double)
           (t (muffle-optimization-notes
                (lla-warn-suboptimal "non-specialized vector")
                ,(funcall body-generator nil))))))))


(defmacro with-vector-type-declarations ((vector &key
                                                 other-vectors
                                                 (simple-test t)
                                                 (vector? t)
                                                 (lla-type (gensym "LLA-TYPE")))
                                         &body body)
  "Like WITH-VECTOR-TYPE-EXPANSION, but constant body-generator taken
from BODY."
  `(with-vector-type-expansion (,vector
                                :other-vectors ,other-vectors
                                :simple-test ,simple-test
                                :vector? ,vector?
                                :lla-type ,lla-type)
     (lambda (lla-type)
       (declare (ignore lla-type))
       `(progn ,',@body))))

(defun common-lla-type (lla-type1 lla-type2)
  "Return a common LLA type that both of the given types can be
coerced to, or NIL if that is not possible.  Arguments have to be
valid LLA types, which is not checked."
  (cond
    ((or (null lla-type1) (null lla-type2))
     nil)
    ((eq lla-type1 lla-type2)
     lla-type1)
    (t (let ((complex? (or (lla-complex? lla-type1)
                           (lla-complex? lla-type2)))
             (double? (or (lla-double? lla-type1)
                          (lla-double? lla-type2))))
         (if complex?
             (if double? :complex-double :complex-single)
             (if double? :double :single))))))

(defgeneric pack (object)
  (:documentation "The purpose of this is to return objects for immediate
  detection of element types, affecting speed favorably.  Vectors are converted
  to simple-arra1.  The returned object may share structure with the
  original."))

(defmethod pack ((vector vector))
  #+(or lla-vector-integer lla-vector-single lla-vector-double
        lla-vector-complex-single lla-vector-complex-double )
  (let ((lla-type (array-lla-type vector)))
    (if lla-type
        (as-simple-array1 vector)
        (let ((common-lla-type 
               (reduce #'common-lla-type vector
                       :key #'atom-representable-lla-type)))
          (if common-lla-type
              (copy-vector vector common-lla-type)
              vector))))
  #+(not (or lla-vector-integer lla-vector-single lla-vector-double
             lla-vector-complex-single lla-vector-complex-double))
  ;; CLISP & ilk: don't bother
  vector)

(define-modify-macro packf () pack "PACK place.")
