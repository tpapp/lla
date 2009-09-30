(in-package :lla)

;;;; ** General notes.
;;;; 
;;;; Currently, everthing is very SBCL-specific, and will remain so
;;;; until things finalize.  In theory, it should be possible to do
;;;; everything in other implementations, albeit more slowly (see the
;;;; FFA for examples).  I (Tamas) decided not to bother about these
;;;; things for a while, concentrating on the API.  This file should
;;;; provide an interface which should not need to be changed for
;;;; other implementations.
;;;;
;;;; What about Rif-style foreign numeric vectors?  Frankly, I see no
;;;; point: given SBCL's pinning ability, foreign numeric vectors do
;;;; not confer any advantage, but have a huge disadvantage: they
;;;; cannot be moved by the GC, and finalization is rather ad-hoc, so
;;;; I see no reason to support them for my purposes (but I am open to
;;;; arguments).

#-sbcl (error "You are not using SBCL!  See note in the source ~
  above.")

;;;; ** numeric vectors
;;;;

(define-abstract-class numeric-vector ()
  ((data :type (simple-array * (*))
	 :initarg :data :reader numeric-vector-data)))

(eval-when (:compile-toplevel :load-toplevel :execute)
  (defun numeric-vector-class (lla-type)
    "Return the name of the numeric vector class corresponding to the
LLA type."
    (check-type lla-type symbol)
    (make-symbol* 'numeric-vector- lla-type)))

;;;;
;;;; printing and element-type
;;;;

(defgeneric nv-element-type (numeric-vector)
  (:documentation "Return the lla-type of the elements in a
  numeric-vector."))

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

;;;; generic xref interface
;;;;
;;;; !! also check for speed? -- Tamas

(defmethod xtype ((nv numeric-vector))
  (lla-type->lisp-type (nv-element-type nv)))

(defmethod xrank ((nv numeric-vector))
  (declare (ignore nv))
  1)

(defmethod xsize ((nv numeric-vector))
  (length (numeric-vector-data nv)))

(defmethod xdim ((nv numeric-vector) axis-number)
  (unless (zerop axis-number)
    (error "numeric-vectors are of rank 1"))
  (xsize nv))

(defmethod xdims ((nv numeric-vector))
  (list (xsize nv)))

(defmethod xref ((nv numeric-vector) &rest subscripts)
  (aref (numeric-vector-data nv) (first subscripts)))

(defmethod (setf xref) (value (nv numeric-vector) &rest subscripts)
  (setf (aref (numeric-vector-data nv) (first subscripts))
	value))

;;;; 
;;;; Define the subclasses.
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
     (defmethod nv-element-type ((numeric-vector ,class-name))
       ',lla-type)
     (defmethod take ((vector vector) (class (eql ',class-name)) &key function 
                      force-copy-p &allow-other-keys)
       (make-instance ',class-name :data
                      (if (and (eq (array-element-type vector)
                                   (upgraded-array-element-type ',lisp-type))
                               (not function)
                               (not force-copy-p))
                          vector
                          (map ',array-type 
                               (if function
                                   (lambda (x) (coerce (funcall function x) ',lisp-type))
                                   (lambda (x) (coerce x ',lisp-type)))
                               vector)))))))

(def-numeric-vector-class :single)
(def-numeric-vector-class :double)
(def-numeric-vector-class :complex-single)
(def-numeric-vector-class :complex-double)
(def-numeric-vector-class :integer)

(defun make-numeric-vector (length lla-type &optional initial-contents)
  "Return a numeric vector with given length and LLA-type.
Initial-contents is interpeted as follows:
  - a number is used to fill the vector, coerced to the correct type
  - a list is copied, elements coerced
  - a vector is coerced then saved, thus it it not copied if the
    elements are of the correct type.
If initial-contents is nil, the array elements are not initialized."
  (check-type lla-type lla-type)
  (check-type length integer)
  (let* ((lisp-type (lla-type->lisp-type lla-type))
	 (data (ctypecase initial-contents
		 (null (make-array length :element-type lisp-type))
		 (number
		    (make-array length :element-type lisp-type
				:initial-element
				(coerce initial-contents lisp-type)))
		 (vector
		    (coerce initial-contents
			    `(simple-array ,lisp-type (,length))))
		 (list
		    (make-array  length :element-type lisp-type
				 :initial-contents
				 (mapcar (lambda (x)
					   (coerce x lisp-type))
					 initial-contents))))))
    (make-instance (numeric-vector-class lla-type)
		   :data data)))

;;;;
;;;;  copying and conversion (between float types)
;;;;

(defgeneric nv-copy (nv type)
  (:documentation "Make a copy of the numeric vector, converting (if
  necessary) to the target type.  Types are LLA type.")
  (:method (nv type)
    ;; fallback method for impossible conversions, signal error
    (error "can't coerce numeric-vector of type ~A to ~A"
	   (nv-element-type nv) type)))

(defmacro def-nv-copy (source-type target-type)
  "Define a method for nv-copy for given source and target
types."
  ;; !!! optimize, mainly with declarations, etc
  (let ((source-lisp-type (lla-type->lisp-type source-type))
	(target-lisp-type (lla-type->lisp-type target-type)))
    `(defmethod nv-copy ((nv ,(numeric-vector-class source-type))
			 (target-type (eql ',target-type)))
       (let* ((data (numeric-vector-data nv))
	      (length (length data))
	      (copy (make-array length :element-type
				',target-lisp-type)))
	 (dotimes (i length)
	   (setf (aref copy i)
		 ,(if (equal source-lisp-type target-lisp-type)
		      '(aref data i)
		      `(coerce (aref data i)
			       ',target-lisp-type))))
	 (make-instance ',(numeric-vector-class target-type)
			:data copy)))))

(for-coercible-pairs (source target)
		     (def-nv-copy source target))

;;;;
;;;;  Implementation of wrapper macros.
;;;;

#+sbcl
(defmacro with-pinned-vector ((vector pointer) &body body)
  "Pin the vector and bind pointer to its data during body."
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
				       (nv-element-type
					,numeric-vector))
				   ,numeric-vector
				   (nv-copy ,numeric-vector
					    ,lla-type))))
	   (with-pinned-vector ((numeric-vector-data ,nv-maybe-copy)
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
	 (let ((,nv-copy (nv-copy ,numeric-vector ,lla-type)))
	   (with-pinned-vector ((numeric-vector-data ,nv-copy)
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
       (let ((,name (make-numeric-vector ,length ,lla-type)))
	 (with-pinned-vector ((numeric-vector-data ,name) ,pointer)
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
       (let ((,name (nv-copy ,numeric-vector ,lla-type)))
	 (with-pinned-vector ((numeric-vector-data ,name) ,pointer)
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
