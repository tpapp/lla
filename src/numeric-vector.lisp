(in-package :lla)

;;;; General notes.

;;; Currently, everthing is very SBCL-specific, and will remain so
;;; until things finalize.  In theory, it should be possible to do
;;; everything in other implementations, albeit more slowly (see the
;;; FFA for examples).  I decided not to bother about these things
;;; for a while, concentrating on the API.  This file should provide
;;; an interface which should not need to be changed for other
;;; implementations. -- Tamas
;;;
;;; What about Rif-style foreign numeric vectors?  Frankly, I see no
;;; point: given SBCL's pinning ability, foreign numeric vectors do
;;; not confer any advantage, but have a disadvantage: they cannot be
;;; moved by the GC, and finalization is rather ad-hoc, so I see no
;;; reason to support them for my purposes (but I am open to
;;; arguments).  Also I rather like the fact that the insides of
;;; NUMERIC-VECTORs can be accessed and manipulated directly as simple
;;; arrays. -- Tamas

;;;; Numeric vectors

;;; 

(define-abstract-class numeric-vector ()
  ((elements :type (simple-array * (*))
	 :initarg :elements :reader elements)
   (shared-p :type boolean :initarg :shared-p :reader shared-p
             :initform nil :documentation "Whether the elements are
shared with another numeric vector."))
  (:documentation "A numeric vector is a wrapper class around a simple
vector ELEMENTS (which may be accessed directly) and SHARED-P, which
is non-nil iff another numeric vector shares the same data.

The semantics of copying, implemented by NV-COPY, is lazy and is
designed to promote functional programming (where possible)."))

;;;; helper functions

(eval-when (:compile-toplevel :load-toplevel :execute)
  (defun numeric-vector-class (lla-type)
    "Return the class name corresponding to LLA-TYPE."
    (append-lla-type numeric-vector lla-type))
  (defun nv-array-type (lla-type &optional length)
  "Return Lisp array type for a NUMERIC-VECTOR of LLA-TYPE."
  `(simple-array ,(upgraded-array-element-type (lla-type->lisp-type lla-type))
                 (,(if length length '*)))))


;;;; some LLA-specific generic interface and utility functions

(defgeneric copy-elements (numeric-vector)
  (:documentation "Replace the elements vector of numeric-vector with
  a copy if it is shared.  Also return the resulting vector.")
  (:method ((nv numeric-vector))
    (with-slots (elements shared-p) nv
      ;; implementation note: we reply on copy-seq being fast
      (when shared-p
        (setf elements (copy-seq elements)))
      elements)))

(defun nv-copy (nv &optional (shared-p t))
  "Copy a numeric vector.  If shared-p is t, data will be shared and
copied later on demand."
  (let ((copy (make-instance (class-of nv) 
                             :elements (elements nv)
                             :shared-p shared-p)))
    (unless shared-p
      (copy-elements copy))
    copy))


(defun make-nv (length lla-type &optional initial-contents
                use-directly-p)
  "Return a numeric vector with given length and LLA-type.
Initial-contents is interpeted as follows:
  - a number is used to fill the vector, coerced to the correct type
  - a list is copied, elements coerced
  - a vector is copied, elements coerced
If initial-contents is nil, the array elements are not necessarily
initialized (depending on the implementation).  If initial-contents is
a list or a vector, then T is accepted as length.

If use-directly-p, then initial-contents is used as data.  It is
checked for type (also length, if not T).

This is designed to be a \"friendly\" interface that should be able to
use any kind of initial contents if they make sense."
  (check-type lla-type lla-type)
  (let* ((lisp-type (lla-type->lisp-type lla-type))
         (length (cond
                   ((not (eq length t)) ; length given
                    (check-type length dimension)
                    length)
                   (t (check-type initial-contents (or list vector))
                    (length initial-contents))))
         (class (numeric-vector-class lla-type))
         (array-type (nv-array-type lla-type length))
         (data (cond
                 (use-directly-p
                  (if (typep initial-contents array-type)
                      initial-contents
                      (error "INITIAL-CONTENTS is not of the required element type or length.")))
                 ((null initial-contents)
                  (make-array length :element-type lisp-type))
                 ((numberp initial-contents)
                  (make-array length :element-type lisp-type
                              :initial-element
                              (coerce initial-contents lisp-type)))
                 ((or (vectorp initial-contents) (listp initial-contents))
                  (map array-type (lambda (x) (coerce x lisp-type))
                       initial-contents))
                 (t (error "couldn't use ~A to initialize a NUMERIC-VECTOR ~
                       of ~A elements " initial-contents length)))))
    (make-instance class :elements data :shared-p nil)))

;;;; 
;;;; define the subclasses
;;;;


(defmacro define-numeric-vector-class (lla-type)
  "!!! documentation"
  (check-type lla-type symbol)
  (let* ((class-name (numeric-vector-class lla-type))
         (lisp-type (lla-type->lisp-type lla-type))
         (array-type (nv-array-type lla-type)))
    `(progn
       ;; class definition
       (defclass ,class-name (numeric-vector)
         ((data :type ,array-type))
         (:documentation ,(format nil "numeric vector of type ~a"
                                  lla-type)))
       ;; LLA-specific interface
       (defmethod lla-type ((numeric-vector ,class-name))
         ',lla-type)
       ;; XARRAY interface (type-specific)
       (defmethod xelttype ((numeric-vector ,class-name))
         ',lisp-type)
       (defmethod take ((class (eql ',class-name)) (vector vector) &key
                        force-copy-p &allow-other-keys)
         (make-nv t ,lla-type vector
                  (and (not force-copy-p)
                       (typep vector ',array-type))))
       (defmethod xcreate ((class (eql ',class-name)) dimensions &optional options)
         (when options (error "This method does not take any options."))
         (bind (((length) dimensions))
           (make-nv length ,lla-type))))))

(define-numeric-vector-class :single)
(define-numeric-vector-class :double)
(define-numeric-vector-class :complex-single)
(define-numeric-vector-class :complex-double)
(define-numeric-vector-class :integer)

;;;; miscellaneous helper functions

(declaim (inline lb-target-type))
(defun lb-target-type (&rest objects)
  "Find common target type of objects to.  Forces floats, should be
used in LAPACK."
  (binary-code->lla-type
   (reduce #'logior objects :key (lambda (object)
                                   (lla-type->binary-code
                                    (lla-type object))))))

(defgeneric nv-convert (nv lla-type)
  (:documentation "Make a copy of the ELEMENTS of a numeric vector,
  converting (if necessary) to the target type.  Types are LLA type.")
  (:method (nv lla-type)
    ;; fallback method for impossible conversions, signal error
    (error "can't coerce numeric-vector of type ~A to ~A"
	   (lla-type nv) lla-type)))

(eval-when (:compile-toplevel :load-toplevel :execute)
  (defun generate-nv-convert (source-type target-type)
    "Define a method for nv-copy-convert for given source and target
types."
    ;; !!! optimize, mainly with declarations, etc
    (let ((source-lisp-type (lla-type->lisp-type source-type))
          (target-lisp-type (lla-type->lisp-type target-type))
          (class (numeric-vector-class source-type)))
      `(defmethod nv-convert ((nv ,class) (target-type (eql ',target-type)))
         ;;       (declare (optimize speed))
         ,(if (equal source-lisp-type target-lisp-type)
              `(copy-seq (elements nv))
              `(let* ((elements (elements nv))
                      (length (length elements))
                      (copy (make-array length :element-type
                                        ',target-lisp-type)))
                 ;; ??? test if map would be faster
                 (dotimes (i length)
                   (setf (aref copy i)
                         (coerce (aref elements i)
                                 ',target-lisp-type)))
                 copy))))))

(defmacro generate-nv-convert-function ()
  `(progn
     ,@(mapcar (lambda (pair)
                 (apply #'generate-nv-convert pair))
               (coercible-pairs-list))))
(generate-nv-convert-function)


;;;; XARRAY interface for NUMERIC-VECTOR

(defgeneric default-lla-type (object)
  (:documentation "Try to determine default lla-type for object.  If
  this cannot be done, return nil.")
  (:method ((vector vector))
    (bind ((vector-lla-type (lisp-type->lla-type 
                             (array-element-type vector)
                             nil)))
      (if vector-lla-type vector-lla-type
          (find-element-type vector))))
  (:method ((array array))
    (default-lla-type (flatten-array array)))
  (:method ((view view))
    (default-lla-type (original-ancestor view)))
  (:method ((nv numeric-vector))
    (lla-type nv)))

(defmethod xcreate ((class (eql 'numeric-vector)) dimensions &optional
                    options)
  (bind (((&key (lla-type :double)) options)
         ((n) dimensions))
    (make-nv n lla-type 0)))

(defmethod take ((class-name (eql 'numeric-vector)) (vector vector) &key force-copy-p
                 options)
  ;; this method tries to guess the type from the provided vector, or
  ;; use a default
  (bind (((&key (lla-type (default-lla-type vector))) options))
    (unless lla-type
      (error "could not determine lla-type"))
    (make-nv t lla-type vector 
             (and (eq lla-type (lisp-type->lla-type
                                (array-element-type vector)))
                  (not force-copy-p)))))

(defmethod take ((class-name (eql 'numeric-vector)) (view view) &key
                 options force-copy-p)
  (declare (ignore force-copy-p))
  (bind (((&key (lla-type (default-lla-type view))) options))
    (unless lla-type
      (error "could not determined lla-type"))
    (bind ((vector (take 'array view
                         :options `(:element-type
                                    ,(lla-type->lisp-type lla-type))
                         :force-copy-p t)))
      (make-nv t lla-type vector t))))

(defmethod take ((class-name (eql 'numeric-vector))
                 (nv numeric-vector) &key options force-copy-p)
  (bind (((&key (lla-type :double)) options))
    (if (eq lla-type (lla-type nv))
        (nv-copy nv (not force-copy-p))
        (nv-convert nv lla-type))))


;; !! also check for speed? -- Tamas

(defmethod xrank ((nv numeric-vector))
  (declare (ignore nv))
  1)

(defmethod xsize ((nv numeric-vector))
  (length (elements nv)))

(defmethod xdim ((nv numeric-vector) axis-number)
  (unless (zerop axis-number)
    (error "numeric-vectors are of rank 1"))
  (xsize nv))

(defmethod xdims ((nv numeric-vector))
  (list (xsize nv)))

(defmethod xref ((nv numeric-vector) &rest subscripts)
  (aref (elements nv) (first subscripts)))

(defmethod (setf xref) (value (nv numeric-vector) &rest subscripts)
  (when (shared-p nv)
    (copy-elements nv))
  (setf (aref (elements nv) (first subscripts))
	value))

;;;; some operations

(defmacro define-nv-elementwise-operation (op &optional (method-name (make-symbol* 'X op)))
  "Define an elementwise operation between vectors name METHOD-NAME."
  `(progn
     ,@(mapcar (lambda (pair)
                 (bind (((a-type b-type) pair)
                        (a-class (numeric-vector-class a-type))
                        (b-class (numeric-vector-class b-type))
                        (a-array-type (nv-array-type a-type))
                        (b-array-type (nv-array-type b-type))
                        (common-type (common-target-type a-type b-type))
                        (common-array-type (nv-array-type common-type)))
                   `(defmethod ,method-name ((a ,a-class) (b ,b-class) &key &allow-other-keys)
                      (let ((a-elements (elements a))
                            (b-elements (elements b)))
                        (declare (type ,a-array-type a-elements)
                                 (type ,b-array-type b-elements))
                        (make-instance ',(numeric-vector-class common-type)
                                       :elements (map ',common-array-type
                                                      (lambda (a-elt b-elt)
                                                        (,op a-elt b-elt))
                                                      a-elements b-elements))))))
               (coercible-pairs-list))))

(define-nv-elementwise-operation +)
(define-nv-elementwise-operation -)
(define-nv-elementwise-operation *)
(define-nv-elementwise-operation /)
