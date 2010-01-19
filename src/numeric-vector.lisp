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
             :initarg :elements :reader elements
             :documentation "Elements, a specialized simple-array.
             Initialized with zeros by default."))
  (:documentation "A numeric vector is a wrapper class around a simple
vector ELEMENTS"))

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

(expand-for-lla-types (lla-type)
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
(defun make-nv* (lla-type elements)
  "Create a NUMERIC-VECTOR of LLA-TYPE, with given ELEMENTS.  Note
that there is no type checking, and elements are not copied: this is
effectively shorthand for a MAKE-INSTANCE call.  For internal use, not
exported."
  (make-instance (nv-class lla-type) :elements elements))

(declaim (inline make-nv-elements))
(defun make-nv-elements (length lla-type &optional (initial-element 0))
  (expand-for-lla-types (lla-type :prologue (ecase lla-type))
    (let ((lisp-type (lla-type->lisp-type lla-type)))
      `(,lla-type (make-array length :element-type ',lisp-type
                              :initial-element (coerce initial-element ',lisp-type))))))
  
  

(defun make-nv (length lla-type &optional (initial-element 0))
  "Create a NUMERIC-VECTOR of LLA-TYPE, optionally initialized with INITIAL-ELEMENTs."
  (make-nv* lla-type (make-nv-elements length lla-type initial-element)))

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

(defun copy-elements-into (source source-type source-index destination destination-type destination-index length)
  "Copy LENGTH elements from SOURCE to DESTINATION, starting at the
given indexes.  Caller promises that SOURCE and DESTINATION are simple
vectors conforming to the corresponding LLA type.  Return no value.

Usage note: call this function whenever you need to copy/convert
elements.  It is supposed to contain the optimized versions for
conversion etc, so nothing else should be optimized."
  ;; !!! need to optimize this sometime
  (declare (ignore source-type))
  (let ((type (lla-type->lisp-type destination-type)))
    (iter
      (for source-i :from source-index :below (+ source-index length))
      (for destination-i :from destination-index)
      (setf (aref destination destination-i) (coerce (aref source source-i) type))))
  (values))

(defun copy-elements (nv &optional (destination-type (lla-type nv)))
  "Return a vector that is a copy if ELEMENTS in NUMERIC-VECTOR,
converting if necessary.  Note: to copy a numeric-vector, just use
COPY-NV."
  (let* ((source (elements nv))
         (length (length source))
         (destination (make-nv-elements length destination-type)))
    (copy-elements-into source (lla-type nv) 0 destination destination-type 0 length)
    destination))

(defun copy-elements% (nv destination-type copy-p)
  "Copy and return elements if source and destination types don't
match, or if COPY-P; otherwise just return elements.  Usage note:
meant to be used in functions that implement the DESTINATION-TYPE &
COPY-P semantics."
  (let ((source-type (lla-type nv)))
    (if (or copy-p (not (eq destination-type source-type)))
        (copy-elements nv destination-type)
        (elements nv))))

(defun float-elements% (nv &optional (float-type :double))
  "If NV has integer elements, convert to FLOAT-TYPE and return,
otherwise just return elements.  Return element type as a second
value.  Usage note: for use in operations which are not closed on
integers (eg /).  *Not exported*."
  (let ((type (lla-type nv)))
    (if (eq type :integer)
        (values (copy-elements nv float-type) float-type)
        (values (elements nv) type))))

(defun copy-nv (nv &key (destination-type (lla-type nv)) (copy-p nil))
  "Copy (or convert) a numeric vector.  If DESTINATION-TYPE is the
same and COPY-P is nil, then ELEMENTS are shared.  If NV is an
instance of a subclass of NUMERIC-VECTOR, the result is still a
NUMERIC-VECTOR."
  (make-instance (nv-class destination-type)
                 :elements (copy-elements% nv destination-type copy-p)))
  
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

(defmethod xcreate ((class (eql 'numeric-vector)) dimensions &optional
                    options)
  (bind (((&key (lla-type :double)) options)
         ((n) dimensions))
    (make-nv n lla-type)))

(expand-for-lla-types (lla-type)
  `(defmethod xcreate ((class (eql ',(nv-class lla-type))) dimensions &optional options)
     (bind (((&key (lla-type ,lla-type)) options)
            ((n) dimensions))
       (assert (eq lla-type ,lla-type) () "Incompatible options.")
       (make-nv n ,lla-type))))

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

(expand-for-lla-types (lla-type)
  `(defmethod xref ((nv ,(nv-class lla-type))  &rest subscripts)
     (declare (optimize (speed 3))
              #+sbcl (sb-ext:muffle-conditions sb-ext:compiler-note))
     (aref (the ,(nv-array-type lla-type) (elements nv)) (first subscripts))))

(expand-for-lla-types (lla-type)
  `(defmethod (setf xref) (value (nv ,(nv-class lla-type)) &rest subscripts)
     (declare (optimize (speed 3))
              #+sbcl (sb-ext:muffle-conditions sb-ext:compiler-note))
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
                                              (result (make-nv-elements length
                                                                        ,common-type)))
                                         (assert (= length (length b-elements)))
                                         (dotimes (i length)
                                           (setf (aref result i)
                         ;;; here the result type is deliberately not
                         ;;; declared, I want to catch overflows
                                                 (,op (aref a-elements i)
                                                      (aref b-elements i))))
                                         result))))))
               (coercible-pairs-list))))

(defmacro define-nv-scalar-operation (op &optional (method-name (make-symbol* 'X op)))
  "Define vector-scalar operation."
  ;;; !!! not optimized, maybe one of these days
  ;;; !!! obviously DOESN'T give the correct type for integers and /
  `(progn
     (defmethod ,method-name ((a numeric-vector) (b number) &key element-type)
       (when element-type
         (error "Not supported."))
       (let* ((elements (elements a))
              (length (length elements))
              (common-type (common-target-type (lla-type a)
                                               (lisp-type->lla-type (type-of b))))
              (result (make-nv-elements length common-type)))
         (dotimes (i length)
           (setf (aref result i) (,op (aref elements i) b)))
         (make-nv* common-type result)))
     (defmethod ,method-name ((a number) (b numeric-vector) &key element-type)
       (when element-type
         (error "Not supported."))
       (let* ((elements (elements b))
              (length (length elements))
              (common-type (common-target-type (lla-type a)
                                               (lisp-type->lla-type (type-of a))))
              (result (make-nv-elements length common-type)))
         (dotimes (i length)
           (setf (aref result i) (,op a (aref elements i))))
         (make-nv* common-type result)))))

(define-nv-elementwise-operation +)
(define-nv-elementwise-operation -)
(define-nv-elementwise-operation *)
(define-nv-elementwise-operation /)

(define-nv-scalar-operation +)
(define-nv-scalar-operation *)
(define-nv-scalar-operation -)
(define-nv-scalar-operation /)
