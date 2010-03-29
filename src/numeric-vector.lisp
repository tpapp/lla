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

;;; Naming convention: class names have NUMERIC-VECTOR spelled out,
;;; but functions/macros abbreviate it as nv.

;;;  Abstract class for numeric vector-like classes (basically all
;;;  classes in LLA).

(define-abstract-class numeric-vector-like ()
  ((lla-type :accessor lla-type :initarg :lla-type)
   (elements :type (simple-array * (*))
             :initarg :elements :reader elements
             :documentation "Elements, a specialized simple-array.
             Initialized with zeros by default."))
  (:documentation "A numeric vector-like object is a wrapper class
around a simple vector ELEMENTS.  This class does not imply that the
elements represent a vector, this is only used as a basis for
objects which can be mapped to a vector somehow."))

(defmethod xelttype ((nv numeric-vector-like))
  (lla-type->lisp-type (lla-type nv)))

(eval-when (:compile-toplevel :load-toplevel :execute)
  (defun nv-array-type (&optional lla-type length)
    "Return Lisp array type for a NUMERIC-VECTOR of LLA-TYPE.  If the
latter is not given, simply return SIMPLE-ARRAY."
    (if lla-type
        `(simple-array ,(upgraded-array-element-type (lla-type->lisp-type lla-type))
                       (,(if length length '*)))
        'simple-array)))

;;;  Numeric vectors

(defclass numeric-vector (numeric-vector-like)
  ()
  (:documentation "Numeric vector, length is the same as of
  ELEMENTS."))

;;;; object creation functions

(declaim (inline make-nv*))
(defun make-nv* (lla-type elements)
  "Create a NUMERIC-VECTOR of LLA-TYPE, with given ELEMENTS.  Note
that there is no type checking, and elements are not copied: this is
effectively shorthand for a MAKE-INSTANCE call.  For internal use, not
exported."
  (make-instance 'numeric-vector :lla-type lla-type :elements elements))

(declaim (inline make-nv-elements))
(defun make-nv-elements (lla-type length &optional (initial-element 0))
  (expand-for-lla-types (lla-type :prologue (ecase lla-type))
    (let ((lisp-type (lla-type->lisp-type lla-type)))
      `(,lla-type (make-array length :element-type ',lisp-type
                              :initial-element (coerce initial-element ',lisp-type))))))

(defun make-nv (lla-type length &optional (initial-element 0))
  "Create a NUMERIC-VECTOR of LLA-TYPE, optionally initialized with INITIAL-ELEMENTs."
  (make-nv* lla-type (make-nv-elements lla-type length initial-element)))

(defun create-nv (initial-contents &optional lla-type)
  "Create numeric vector with given initial contents (a list or a
vector).  Unless LLA-TYPE is given, it is inferred from the elements.
This is a convenience function for easily creation of NUMERIC-VECTORs.
Also see *FORCE-FLOAT* and *FORCE-DOUBLE*."
  (let* ((lla-type (infer-lla-type lla-type initial-contents))
         (lisp-type (lla-type->lisp-type lla-type))
         (length (length initial-contents)))
    (make-nv* lla-type (map (nv-array-type lla-type length)
                            (lambda (x) (coerce x lisp-type))
                            initial-contents))))

;;;; load forms

(defun make-elements-load-form (nv)
  (bind (((:slots-read-only elements) nv))
    `(make-array ,(length elements)
                 :element-type ',(lla-type->lisp-type (lla-type nv))
                 :initial-contents ,elements)))

(defmethod make-load-form ((nv numeric-vector) &optional environment)
  (declare (ignore environment))
  `(make-nv* ,(lla-type nv) ,(make-elements-load-form nv)))

;;;; some LLA-specific generic interface and utility functions

(defmacro copy-elements-loop% (&optional (source-element-form 
                                          '(row-major-aref source source-index)))
  "Loop for copy copying elements.  Used in copy-elements etc, see
that function for variable names."
  `(iter
     (for source-index :from source-offset :below (+ source-offset length))
     (for destination-index :from destination-offset)
     (declare (iterate:declare-variables)
              (fixnum source-index destination-index))
     (setf (row-major-aref destination destination-index)
           ,source-element-form)))

(defun copy-elements% (length source source-offset source-type
                       destination destination-offset)
  "Fast version of copy-elements for coinciding types."
  (declare (optimize (speed 3) (safety 0))
           (fixnum length source-offset destination-offset))
  (expand-for-lla-types (lla-type :prologue (ecase source-type))
    `(,lla-type
      (locally
          (declare (type ,(nv-array-type lla-type) source destination))
        (copy-elements-loop%)))))

(defun copy-elements (length source source-offset source-type
                      destination destination-offset
                      &optional (destination-type source-type))
  "Copy LENGTH elements from SOURCE to DESTINATION, starting at the
given offsets.  Caller promises that SOURCE and DESTINATION are simple
vectors conforming to the corresponding LLA type.  Return no value.

Usage note: call this function whenever you need to copy/convert
elements.  It is supposed to contain the optimized versions for
conversion etc, so nothing else should be optimized."
  (if (eq source-type destination-type)
      (copy-elements% length source source-offset source-type
                      destination destination-offset)
      (let ((destination-lisp-type (lla-type->lisp-type destination-type)))
        (copy-elements-loop% (coerce (row-major-aref source source-index)
                                    destination-lisp-type))))
  (values))

(defun copy-nv-elements (nv &key (destination-type (lla-type nv))
                      (length (length (elements nv))))
  "Return a vector that is a copy if ELEMENTS in NUMERIC-VECTOR,
converting if necessary.  Note: to copy a numeric-vector, just use
COPY-NV."
  (let* ((source (elements nv))
         (destination (make-nv-elements destination-type length)))
    (copy-elements length source 0 (lla-type nv) 
                   destination 0 destination-type)
    destination))

(defun copy-nv-elements% (nv destination-type copy-p)
  "Copy and return elements if source and destination types don't
match, or if COPY-P; otherwise just return elements.  Usage note:
meant to be used in functions that implement the DESTINATION-TYPE &
COPY-P semantics."
  (let ((source-type (lla-type nv)))
    (if (or copy-p (not (eq destination-type source-type)))
        (copy-nv-elements nv :destination-type destination-type)
        (elements nv))))

(defun float-elements% (nv &optional (float-type :double))
  "If NV has integer elements, convert to FLOAT-TYPE and return,
otherwise just return elements.  Return element type as a second
value.  Usage note: for use in operations which are not closed on
integers (eg /).  *Not exported*."
  (let ((type (lla-type nv)))
    (if (eq type :integer)
        (values (copy-nv-elements nv :destination-type float-type) float-type)
        (values (elements nv) type))))

(defun copy-nv (nv &key (destination-type (lla-type nv)) (copy-p nil))
  "Copy (or convert) a numeric vector.  If DESTINATION-TYPE is the
same and COPY-P is nil, then ELEMENTS are shared.  If NV is an
instance of a subclass of NUMERIC-VECTOR, the result is still a
NUMERIC-VECTOR."
  (make-instance 'numeric-vector :lla-type (lla-type nv)
                 :elements (copy-nv-elements% nv destination-type copy-p)))
  
;;; XARRAY interface

(defmethod xrank ((nv numeric-vector-like))
  (declare (ignore nv))
  1)

(defmethod xsize ((nv numeric-vector-like))
  (length (elements nv)))

(defmethod xdim ((nv numeric-vector-like) axis-number)
  (assert (zerop axis-number) () 'xdim-invalid-axis-number)
  (xsize nv))

(defmethod xdims ((nv numeric-vector-like))
  (list (xsize nv)))

(defmethod xref ((nv numeric-vector-like) &rest subscripts)
  (assert (not (cdr subscripts)) () 'xref-wrong-number-of-subscripts)
  (aref (elements nv) (first subscripts)))

(defmethod (setf xref) (value (nv numeric-vector-like) &rest subscripts)
  (assert (not (cdr subscripts)) () 'xref-wrong-number-of-subscripts)
  (setf (aref (elements nv) (first subscripts))
        value))
