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

(defclass numeric-vector ()
  ((lla-type :accessor lla-type :initarg :lla-type)
   (elements :type (simple-array * (*))
             :initarg :elements :reader elements
             :documentation "Elements, a specialized simple-array.
             Initialized with zeros by default."))
  (:documentation "A numeric vector is a wrapper class around a simple
vector ELEMENTS"))

;;;; object creation functions

(defun nv-array-type (&optional lla-type length)
  "Return Lisp array type for a NUMERIC-VECTOR of LLA-TYPE.  If the
latter is not given, simply return SIMPLE-ARRAY."
  (if lla-type
      `(simple-array ,(upgraded-array-element-type (lla-type->lisp-type lla-type))
                     (,(if length length '*)))
      'simple-array))

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

(defun copy-elements-into (source source-type source-index
                           destination destination-type destination-index
                           length)
  "Copy LENGTH elements from SOURCE to DESTINATION, starting at the
given indexes.  Caller promises that SOURCE and DESTINATION are simple
vectors conforming to the corresponding LLA type.  Return no value.

Usage note: call this function whenever you need to copy/convert
elements.  It is supposed to contain the optimized versions for
conversion etc, so nothing else should be optimized."
  ;; !!! need to optimize this some day
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
         (destination (make-nv-elements destination-type length)))
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
  (make-instance 'numeric-vector :lla-type (lla-type nv)
                 :elements (copy-elements% nv destination-type copy-p)))
  
;;; XARRAY interface

(defmethod xelttype ((nv numeric-vector))
  (lla-type->lisp-type (lla-type nv)))

(defmethod xrank ((nv numeric-vector))
  (declare (ignore nv))
  1)

(defmethod xsize ((nv numeric-vector))
  (length (elements nv)))

(defmethod xdim ((nv numeric-vector) axis-number)
  (assert (zerop axis-number) () 'xdim-invalid-axis-number)
  (xsize nv))

(defmethod xdims ((nv numeric-vector))
  (list (xsize nv)))

(defmethod xref ((nv numeric-vector) &rest subscripts)
  (assert (not (cdr subscripts)) () 'xref-wrong-number-of-subscripts)
  (aref (elements nv) (first subscripts)))

(defmethod (setf xref) (value (nv numeric-vector) &rest subscripts)
  (assert (not (cdr subscripts)) () 'xref-wrong-number-of-subscripts)
  (setf (aref (elements nv) (first subscripts))
        value))
