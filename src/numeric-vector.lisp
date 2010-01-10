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

;;; Naming convention: class names have NUMERIC-VECTOR spelled out,
;;; but functions/macros abbreviate it as nv.

(define-abstract-class numeric-vector ()
  ((elements :type (simple-array * (*))
	 :initarg :elements :reader elements)
   (shared-p :type boolean :initarg :shared-p :reader shared-p
             :initform nil :documentation "Whether the elements are
shared with another numeric vector."))
  (:documentation "A numeric vector is a wrapper class around a simple
vector ELEMENTS (which may be accessed directly) and SHARED-P, which
is non-nil iff another numeric vector shares the same data.

The semantics of copying, implemented by COPY-NV, is lazy and is
designed to promote functional programming (where possible)."))

;;;; helper functions

(eval-when (:compile-toplevel :load-toplevel :execute)
  (defun nv-class (lla-type)
    "Return the class name corresponding to LLA-TYPE."
    (append-lla-type numeric-vector lla-type))
  (defun nv-array-type (lla-type &optional length)
  "Return Lisp array type for a NUMERIC-VECTOR of LLA-TYPE."
  `(simple-array ,(upgraded-array-element-type (lla-type->lisp-type lla-type))
                 (,(if length length '*)))))


;;;; define the subclasses

(expand-for-lla-types lla-type
  (let* ((class-name (nv-class lla-type))
         (lisp-type (lla-type->lisp-type lla-type))
         (array-type (nv-array-type lla-type)))
    `(progn
       ;; class definition
       (defclass ,class-name (numeric-vector)
         ((data :type ,array-type))
         (:documentation ,(format nil "numeric vector of type ~a"
                                  lla-type)))
       ;; LLA-type
       (defmethod lla-type ((numeric-vector ,class-name))
         ',lla-type)
       ;; XELTTYPE
       (defmethod xelttype ((numeric-vector ,class-name))
         ',lisp-type))))


;;;; some LLA-specific generic interface and utility functions

(declaim (inline make-nv*))
(defun make-nv* (lla-type elements &optional (shared-p nil))
  "Create a NUMERIC-VECTOR of LLA-TYPE, with given ELEMENTS.  Note
that there is no type checking, and elements are not copied: this is
effectively shorthand for a MAKE-INSTANCE call.  For internal use, not
exported."
  (make-instance (nv-class lla-type) :elements elements :shared-p shared-p))

(defun make-nv (length lla-type &optional initial-element)
  "Create a NUMERIC-VECTOR of LLA-TYPE, optionally initialized with INITIAL-ELEMENTs."
  (let ((lisp-type (lla-type->lisp-type lla-type)))
    (make-nv* lla-type (if initial-element
                           (make-array length :element-type lisp-type
                                       :initial-element (coerce initial-element lisp-type))
                           (make-array length :element-type lisp-type)))))

(defun create-nv (initial-contents &optional lla-type)
  "Create numeric vector with given initial contents (a list or a
vector).  Unless LLA-TYPE is given, it is inferred from the elements.
This is a convenience function for easily creation of NUMERIC-VECTORs.
Also see *forced-float*."
  (let* ((lla-type (infer-lla-type lla-type initial-contents *force-float*))
         (lisp-type (lla-type->lisp-type lla-type))
         (length (length initial-contents)))
    (make-nv* lla-type (map (nv-array-type lla-type length)
                            (lambda (x) (coerce x lisp-type))
                            initial-contents))))

(defgeneric copy-elements (nv)
  (:documentation "Return a vector that is a copy if ELEMENTS in
NUMERIC-VECTOR.  Note: to copy a numeric-vector, just use COPY-NV,
which can make a lazy copy."))
(expand-for-lla-types lla-type
  `(defmethod copy-elements ((nv ,(nv-class lla-type)))
     (declare (optimize speed))
     (bind (((:slots-read-only elements) nv))
       (declare (type ,(nv-array-type lla-type) elements))
       (copy-seq elements))))

(defun ensure-unshared (nv)
  "Replace the elements vector of numeric-vector with a copy if it is
shared.  Return no values."
    (with-slots (elements shared-p) nv
      ;; implementation note: we rely on copy-seq being fast
      (when shared-p
        (setf shared-p nil)
        (setf elements (copy-elements nv)))
      (values)))

(defgeneric copy-nv (nv)
  (:documentation "Copy a numeric vector.  ELEMENTS are shared.  If NV
is an instance of a subclass of NUMERIC-VECTOR, the result is still a
NUMERIC-VECTOR."))
(expand-for-lla-types lla-type
  (let ((class (nv-class lla-type)))
    `(defmethod copy-nv ((nv ,class))
       (make-instance ',class
                      :elements (elements nv)))))

(defgeneric convert-elements (nv lla-type)
  (:documentation "Make a (non-shared) copy of the ELEMENTS of a
  numeric vector, converting (if necessary) to the target type.  Types
  are LLA type.")
  (:method (nv lla-type)
    ;; fallback method for impossible conversions, signal error
    (error "can't coerce numeric-vector of type ~A to ~A"
	   (lla-type nv) lla-type)))

(eval-when (:compile-toplevel :load-toplevel :execute)
  (defun generate-convert-elements (source-type target-type)
    "Define a method for convert-nv for given source and target
types."
    (let ((source-lisp-type (lla-type->lisp-type source-type))
          (target-lisp-type (lla-type->lisp-type target-type))
          (class (nv-class source-type)))
      `(defmethod convert-elements ((nv ,class) (target-type (eql ,target-type)))
         (declare (optimize speed)
                  #+sbcl (sb-ext:muffle-conditions sb-ext:compiler-note))
         ,(if (equal source-lisp-type target-lisp-type)
              `(copy-elements nv)
              `(let ((elements (elements nv)))
                 (declare (type ,(nv-array-type source-type) elements))
                 (let* ((length (length elements))
                        (result (make-array length :element-type ',target-lisp-type)))
                   (dotimes (i length)
                     (setf (aref result i)
                           (coerce (aref elements i) ',target-lisp-type)))
                   result)))))))

(defmacro define-convert-elements ()
  `(progn
     ,@(mapcar (lambda (pair)
                 (apply #'generate-convert-elements pair))
               (coercible-pairs-list))))
(define-convert-elements)


;;;; XARRAY interface for NUMERIC-VECTOR

;; (defgeneric default-lla-type (object)
;;   (:documentation "Try to determine default lla-type for object.  If
;;   this cannot be done, return nil.")
;;   (:method ((vector vector))
;;     (bind ((vector-lla-type (lisp-type->lla-type 
;;                              (array-element-type vector)
;;                              nil)))
;;       (if vector-lla-type vector-lla-type
;;           (find-element-type vector))))
;;   (:method ((array array))
;;     (default-lla-type (make-array (array-total-size array)
;;                                   :element-type (array-element-type array)
;;                                   :displaced-to array)))
;;   (:method ((view view))
;;     (default-lla-type (original-ancestor view)))
;;   (:method ((nv numeric-vector))
;;     (lla-type nv)))

;; (defmethod xcreate ((class (eql 'numeric-vector)) dimensions &optional
;;                     options)
;;   (bind (((&key (lla-type :double)) options)
;;          ((n) dimensions))
;;     (make-nv n lla-type)))

;; (defmethod take ((class-name (eql 'numeric-vector)) (vector vector) &key force-copy-p
;;                  options)
;;   ;; this method tries to guess the type from the provided vector, or
;;   ;; use a default
;;   (declare (ignore force-copy-p))
;;   (bind (((&key (lla-type (default-lla-type vector))) options))
;;     (unless lla-type
;;       (error "could not determine lla-type"))
;;     (make-nv nil lla-type vector)))
             

;; (defmethod take ((class-name (eql 'numeric-vector)) (view view) &key
;;                  options force-copy-p)
;;   (declare (ignore force-copy-p))
;;   (bind (((&key (lla-type (default-lla-type view))) options))
;;     (unless lla-type
;;       (error "could not determined lla-type"))
;;     (bind ((vector (take 'array view
;;                          :options `(:element-type
;;                                     ,(lla-type->lisp-type lla-type))
;;                          :force-copy-p t)))
;;       (make-instance (nv-class lla-type) :elements vector))))

;; (defmethod take ((class-name (eql 'numeric-vector))
;;                  (nv numeric-vector) &key options force-copy-p)
;;   (bind (((&key (lla-type :double)) options))
;;     (if (eq lla-type (lla-type nv))
;;         (copy-nv nv (not force-copy-p))
;;         (convert-nv nv lla-type))))


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

(expand-for-lla-types lla-type
  `(defmethod xref ((nv ,(nv-class lla-type))  &rest subscripts)
     (declare (optimize (speed 3))
              #+sbcl (sb-ext:muffle-conditions sb-ext:compiler-note))
     (aref (the ,(nv-array-type lla-type) (elements nv)) (first subscripts))))

(expand-for-lla-types lla-type
  `(defmethod (setf xref) (value (nv ,(nv-class lla-type)) &rest subscripts)
     (declare (optimize (speed 3))
              #+sbcl (sb-ext:muffle-conditions sb-ext:compiler-note))
     (ensure-unshared nv)
     (setf (aref (the ,(nv-array-type lla-type) (elements nv)) (first subscripts))
           value)))

;;;; some operations

(defmacro define-nv-elementwise-operation (op &optional (method-name (make-symbol* 'X op)))
  "Define an elementwise operation between vectors name METHOD-NAME."
  `(progn
     ,@(mapcar (lambda (pair)
                 (bind (((a-type b-type) pair)
                        (a-class (nv-class a-type))
                        (b-class (nv-class b-type))
                        (a-array-type (nv-array-type a-type))
                        (b-array-type (nv-array-type b-type))
                        (common-type (common-target-type a-type b-type)))
                   `(defmethod ,method-name ((a ,a-class) (b ,b-class)
                                             &key &allow-other-keys)
                      (declare (optimize speed)
                               #+sbcl (sb-ext:muffle-conditions sb-ext:compiler-note))
                      (let ((a-elements (elements a))
                            (b-elements (elements b)))
                        (declare (type ,a-array-type a-elements)
                                 (type ,b-array-type b-elements))
                        (make-instance ',(nv-class common-type)
                                       :elements
                                       (let* ((length (length a-elements))
                                              (result (make-array length :element-type 
                                                                  ',(lla-type->lisp-type 
                                                                     common-type))))
                                         (assert (= length (length b-elements)))
                                         (dotimes (i length)
                                           (setf (aref result i)
                         ;;; here the result type is deliberately not
                         ;;; declared, I want to catch overflows
                                                 (,op (aref a-elements i)
                                                      (aref b-elements i))))
                                         result))))))
               (coercible-pairs-list))))

(define-nv-elementwise-operation +)
(define-nv-elementwise-operation -)
(define-nv-elementwise-operation *)
(define-nv-elementwise-operation /)
