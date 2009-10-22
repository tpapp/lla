(in-package :lla)

;;;; ** General notes.
;;;; 
;;;; Currently, everthing is very SBCL-specific, and will remain so
;;;; until things finalize.  In theory, it should be possible to do
;;;; everything in other implementations, albeit more slowly (see the
;;;; FFA for examples).  I decided not to bother about these things
;;;; for a while, concentrating on the API.  This file should provide
;;;; an interface which should not need to be changed for other
;;;; implementations. -- Tamas
;;;;
;;;; What about Rif-style foreign numeric vectors?  Frankly, I see no
;;;; point: given SBCL's pinning ability, foreign numeric vectors do
;;;; not confer any advantage, but have a huge disadvantage: they
;;;; cannot be moved by the GC, and finalization is rather ad-hoc, so
;;;; I see no reason to support them for my purposes (but I am open to
;;;; arguments).  Also I rather like the fact that the insides of
;;;; NUMERIC-VECTORs can be accessed and manipulated directly as
;;;; simple arrays. -- Tamas

;;;; NOTE: everything can and will be made to work with other
;;;; implementations, but I am concentrating on SBCL for the moment.

#-sbcl (error "You are not using SBCL!  See note in the source ~
  above.")

;;;; ** numeric vectors
;;;;
;;;; A numeric vector is a wrapper class around a simple vector DATA
;;;; (which may be accessed directly) and SHARED-P, which is non-nil
;;;; iff another numeric vector shares the same data.
;;;;
;;;; The semantics of copying is lazy and is designed to promote
;;;; functional programming when possible.  nv-copy, which should be
;;;; used for copying NUMERIC-VECTORs

(define-abstract-class numeric-vector ()
  ((data :type (simple-array * (*))
	 :initarg :data :reader nv-data)
   (shared-p :type boolean :initarg :shared-p :reader shared-p
             :initform nil :documentation "Whether the data is
shared with another numeric vector.")))

(eval-when (:compile-toplevel :load-toplevel :execute)
  (defun numeric-vector-class (lla-type)
    "Return the name of the numeric vector class corresponding to the
LLA type."
    (check-type lla-type symbol)
    (make-symbol* 'numeric-vector- lla-type)))

;;;;
;;;; some LLA-specific generic interface and utility functions
;;;;

(defgeneric lla-type (numeric-vector)
  (:documentation "Return the lla-type of the elements of the
  object."))

(defgeneric copy-data (numeric-vector)
  (:documentation "Replace the data vector of numeric-vector with a
  copy if it is shared.  Also return the resulting vector.")
  (:method ((nv numeric-vector))
    (with-slots (data shared-p) nv
      ;; implementation note: we reply on copy-seq being fast
      (when shared-p
        (setf data (copy-seq data)))
      data)))


;;;;
;;;; printing
;;;;

(defmethod print-object ((obj numeric-vector) stream)
  (print-unreadable-object (obj stream :type t)
    (with-slots (data) obj
      (let* ((length (length data))
	     (truncated-length (min *print-length* length)))
	(format stream "of ~a elements: " length)
	(dotimes (i truncated-length)
	  (princ (standard-numeric-formatter (aref data i)) stream)
	  (when (< (1+ i) truncated-length)
	    (princ #\space stream)))
	(when (< truncated-length length)
	  (princ " ..." stream))))))

;;;; 
;;;; define the subclasses
;;;;

(defmacro def-numeric-vector-class (lla-type)
  "!!! documentation"
  (check-type lla-type symbol)
  (let* ((class-name (numeric-vector-class lla-type))
         (lisp-type (lla-type->lisp-type lla-type))
         (array-type `(simple-array ,lisp-type (*))))
  `(progn
     (defclass ,class-name (numeric-vector)
       ((data :type ,array-type))
       (:documentation ,(format nil "numeric vector of type ~a"
				lla-type)))
     (defmethod lla-type ((numeric-vector ,class-name))
       ',lla-type)
     (defmethod xelttype ((numeric-vector ,class-name))
       ',lisp-type)
     (defmethod take ((class (eql ',class-name)) (vector vector) &key function 
                      force-copy-p &allow-other-keys)
       (make-instance ',class-name :data
                      (if (and (eq (array-element-type vector)
                                   (upgraded-array-element-type ',lisp-type))
                               (not function)
                               (not force-copy-p))
                          vector
                          (map ',array-type 
                               (if function
                                   (lambda (x) (coerce (funcall function x)
                                                       ',lisp-type))
                                   (lambda (x) (coerce x ',lisp-type)))
                               vector))))
     (defmethod xcreate ((class (eql ',class-name)) dimensions &key &allow-other-keys)
       (let ((length (cond
                    ((atom dimensions) dimensions)
                    (t (bind (((length &rest rest) dimensions))
                         (when rest
                           (error "~A can only accept a single dimension" ',class-name))
                         length)))))
         (make-nv length ,lla-type))))))

(def-numeric-vector-class :single)
(def-numeric-vector-class :double)
(def-numeric-vector-class :complex-single)
(def-numeric-vector-class :complex-double)
(def-numeric-vector-class :integer)

(defun make-nv (length lla-type &optional initial-contents
                use-directly-p)
  "Return a numeric vector with given length and LLA-type.
Initial-contents is interpeted as follows:
  - a number is used to fill the vector, coerced to the correct type
  - a list is copied, elements coerced
  - a vector is copied, elements coerced
If initial-contents is nil, the array elements are not initialized.

If use-directly-p, then initial-contents is used as data.  It is
checked for type (incl length).

This is designed to be a \"friendly\" interface that should be able to
use any kind of initial contents if they make sense."
  (check-type lla-type lla-type)
  (check-type length integer)
  (let* ((lisp-type (lla-type->lisp-type lla-type))
         (data
          (cond
            (use-directly-p 
             (let ((expected-type `(simple-array
                                    ,(upgraded-array-element-type lisp-type)
                                    (,length))))
               (if (typep initial-contents expected-type)
                   initial-contents
                   (error "INITIAL-CONTENTS ~A is not of the required ~
                           type ~A." initial-contents expected-type))))
            ((null initial-contents)
             (make-array length :element-type lisp-type))
            ((numberp initial-contents)
             (make-array length :element-type lisp-type
                         :initial-element
                         (coerce initial-contents lisp-type)))
            ((or (vectorp initial-contents) (listp initial-contents))
             (map `(simple-array ,lisp-type (,length))
                  (lambda (x) (coerce x lisp-type))
                  initial-contents))
            (t (error "couldn't use ~A to initialize a NUMERIC-VECTOR ~
                       of ~A elements " initial-contents length)))))
    (make-instance (numeric-vector-class lla-type)
		   :data data :shared-p nil)))

;;;; generic xref interface
;;;;
;;;; !! also check for speed? -- Tamas

(defmethod xrank ((nv numeric-vector))
  (declare (ignore nv))
  1)

(defmethod xsize ((nv numeric-vector))
  (length (nv-data nv)))

(defmethod xdim ((nv numeric-vector) axis-number)
  (unless (zerop axis-number)
    (error "numeric-vectors are of rank 1"))
  (xsize nv))

(defmethod xdims ((nv numeric-vector))
  (list (xsize nv)))

(defmethod xref ((nv numeric-vector) &rest subscripts)
  (aref (nv-data nv) (first subscripts)))

(defmethod (setf xref) (value (nv numeric-vector) &rest subscripts)
  (when (shared-p nv)
    (copy-data nv))
  (setf (aref (nv-data nv) (first subscripts))
	value))

;;;;
;;;;  copying and conversion (between float types)
;;;;

(defgeneric nv-copy-convert (nv lla-type)
  (:documentation "Make a copy of the numeric vector, converting (if
  necessary) to the target type.  Types are LLA type.")
  (:method (nv lla-type)
    ;; fallback method for impossible conversions, signal error
    (error "can't coerce numeric-vector of type ~A to ~A"
	   (lla-type nv) lla-type)))

(defmacro def-nv-copy-convert (source-type target-type)
  "Define a method for nv-copy-convert for given source and target
types."
  ;; !!! optimize, mainly with declarations, etc
  (let ((source-lisp-type (lla-type->lisp-type source-type))
	(target-lisp-type (lla-type->lisp-type target-type)))
    `(defmethod nv-copy-convert ((nv ,(numeric-vector-class source-type))
                                 (target-type (eql ',target-type)))
       (let* ((data (nv-data nv))
	      (length (length data))
	      (copy (make-array length :element-type
				',target-lisp-type)))
         ;;; ??? test if map would be faster
	 (dotimes (i length)
	   (setf (aref copy i)
		 ,(if (equal source-lisp-type target-lisp-type)
		      '(aref data i)
		      `(coerce (aref data i)
			       ',target-lisp-type))))
	 (make-instance ',(numeric-vector-class target-type)
			:data copy)))))

(for-coercible-pairs (source target)
		     (def-nv-copy-convert source target))

(defun nv-copy (nv &optional (shared-p t))
  "Copy a numeric vector.  If shared-p is t, data will be shared and
copied later on demand."
  (let ((copy (make-instance (class-of nv) :data (nv-data nv) :shared-p
                             shared-p)))
    (unless shared-p
      (copy-data copy))
    copy))

;;;;
;;;;  Implementation of wrapper macros.
;;;;

#+sbcl
(defmacro with-pinned-vector ((vector pointer) &body body)
  "Pin the vector and bind pointer to its data during body.  This is a
utility function for implementations with pinning, and should not be
used directly in other files, ie it is not part of the interface."
  (check-type pointer symbol)
  (once-only (vector)
    `(sb-sys:with-pinned-objects (,vector)
       (let ((,pointer (sb-sys:vector-sap ,vector)))
	 ,@body))))

#+sbcl
(defmacro with-nv-input ((numeric-vector pointer lla-type) &body body)
  "Makes sure that the contents of numeric-vector are available at
pointer for the duration of body, with the type lla-type (converting
if necessary).  The body should NOT change the data at the pointer in
any way, it is for reading only."
  (check-type pointer symbol)
  (once-only (numeric-vector lla-type)
    (with-unique-names (nv-maybe-copy)
      `(progn
	 (check-type ,numeric-vector numeric-vector)
	 (check-type ,lla-type lla-type)
	 (let ((,nv-maybe-copy (if (eq ,lla-type 
				       (lla-type ,numeric-vector))
				   ,numeric-vector
				   (nv-copy-convert ,numeric-vector
                                                    ,lla-type))))
	   (with-pinned-vector ((nv-data ,nv-maybe-copy)
				,pointer)
	     ,@body))))))

#+sbcl
(defmacro with-nv-input-copied ((numeric-vector pointer lla-type) &body body)
  "Makes sure that the contents of numeric-vector are available at
pointer for the duration of body, with the type lla-type (converting
if necessary).  Data is copied, so it can be changed, but the result
is discarded."
  (check-type pointer symbol)
  (once-only (numeric-vector lla-type)
    (with-unique-names (nv-copy)
      `(progn
	 (check-type ,numeric-vector numeric-vector)
	 (check-type ,lla-type lla-type)
	 (let ((,nv-copy (nv-copy-convert ,numeric-vector ,lla-type)))
	   (with-pinned-vector ((nv-data ,nv-copy)
				,pointer)
	     ,@body))))))

#+sbcl
(defmacro with-nv-output ((name length pointer lla-type) &body body)
  "Allocates a memory area of the given type and length for the
duration of body, and makes sure that the contents are assigned to the
variable named name at the end."
  (check-type name symbol)
  (check-type pointer symbol)
  (once-only (lla-type)
    `(progn
       (check-type ,lla-type lla-type)
       (let ((,name (make-nv ,length ,lla-type)))
	 (with-pinned-vector ((nv-data ,name) ,pointer)
	   ,@body)))))

#+sbcl
(defmacro with-nv-input-output ((numeric-vector name pointer lla-type)
				&body body)
  "A combination of with-nv-input and -output: input is always copied
\(and converted if necessary), is available during body, and the
contents end up in a numeric vector of lla-type, assigned to name."
  (check-type name symbol)
  (check-type pointer symbol)
  (once-only (numeric-vector lla-type)
    `(progn
       (check-type ,numeric-vector numeric-vector)
       (check-type ,lla-type lla-type)
       (let ((,name (nv-copy-convert ,numeric-vector ,lla-type)))
	 (with-pinned-vector ((nv-data ,name) ,pointer)
	   ,@body)))))

(defmacro with-work-area ((pointer lla-type size) &body body)
  "Allocate a work area of size lla-type elements during body,
assigning the pointer to pointer."
  (check-type pointer symbol)
  `(with-foreign-pointer (,pointer 
			  (* ,size (foreign-size* ,lla-type)))
     ,@body))


;;;; 
;;;;  Some LAPACK routines return real and imaginary parts of vectors
;;;;  separately, we have to assemble them.
;;;;

(defun nv-zip-complex-double (pointer n &optional check-real-p)
  "Return the complex numbers stored at pointer (n real parts,
followed by n imaginary parts) as a numeric vector (either double or
complex-double).  If check-real-p, then check if the imaginary part is
0 and if so, return a numeric-vector-double, otherwise always return a
complex-double one."
  (let ((real-p (and check-real-p 
                     (iter
                       (for i :from 0 :below n)
                       (always (zerop (mem-aref pointer :double (+ n i))))))))
    (if real-p
        (let ((data (make-array n :element-type 'double-float)))
          (iter
            (for i :from 0 :below n)
            (setf (aref data i) (mem-aref pointer :double i)))
          (make-instance 'numeric-vector-double :data data))
        (let ((data (make-array n :element-type '(complex double-float))))
          (iter
            (for i :from 0 :below n)
            (setf (aref data i) (complex (mem-aref pointer :double i)
                                         (mem-aref pointer :double (+ n i)))))
          (make-instance 'numeric-vector-complex-double :data data)))))
