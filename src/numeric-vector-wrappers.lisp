(in-package :lla)

;;;; Wrapper macros

;;; Wrapper macros map a numeric vector to a memory area of value wich
;;; specified LLA-TYPE, either for input, input-output or output,
;;; converting the input if required.

;;; NOTE: everything can and will be made to work with other
;;; implementations, but I am concentrating on SBCL for the moment.

#-sbcl (error "You are not using SBCL!  See note in the source ~
  above.")


;;;; All-in-one wrapper macro for numeric vectors.  Just calls one of
;;;; the three with-input-* macros below.  Implementation specific
;;;; code should just modify the three macros below.

(defmacro with-nv-input (((numeric-vector &optional keyword output) pointer lla-type)
                         &body body)
  "Make sure that the contents of NUMERIC-VECTOR are available at
pointer for the duration of BODY.  Optional argument syntax:
  (numeric-vector) is when body promises not to change the elements,
  (numeric-vector :copied) copies the elements into a new location,
  (numeric-vector :output-to output) copies and makes the output
  available as a lisp vector of conforming type."
  (cond
    ((and (symbolp output) (eq keyword :output-to))
     ;; (numeric-vector :output-to output)
     `(with-nv-input-output (,numeric-vector ,output ,pointer ,lla-type)
        ,@body))
    ((and (eq keyword :copied) (null output))
     `(with-nv-input-copied (,numeric-vector ,pointer ,lla-type)
        ,@body))
    ((null keyword)
     `(with-nv-input-readonly (,numeric-vector ,pointer ,lla-type)
        ,@body))
    (t (error "Invalid specification."))))


;;;;
;;;;  Implementation of wrapper macros.
;;;;

;;; The SBCL implementation just uses pinning, and NV-CONVERT for
;;; conversion.

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
(defmacro with-nv-input-readonly ((numeric-vector pointer lla-type) &body body)
  "Makes sure that the contents of numeric-vector are available at
pointer for the duration of body, with the type lla-type (converting
if necessary).  The body should NOT change the data at the pointer in
any way, it is for reading only."
  (check-type pointer symbol)
  (once-only (numeric-vector lla-type)
    (with-unique-names (elements)
      `(progn
	 (check-type ,numeric-vector numeric-vector)
	 (check-type ,lla-type lla-type)
	 (let ((,elements (if (eq ,lla-type 
                                  (lla-type ,numeric-vector))
                              (copy-elements ,numeric-vector)
                              (convert-elements ,numeric-vector ,lla-type))))
	   (with-pinned-vector (,elements ,pointer)
	     ,@body))))))

#+sbcl
(defmacro with-nv-input-copied ((numeric-vector pointer lla-type) &body body)
  "Makes sure that the contents of numeric-vector are available at
pointer for the duration of body, with the type lla-type (converting
if necessary).  Data is copied, so it can be changed, but the result
is discarded."
  (check-type pointer symbol)
  (once-only (numeric-vector lla-type)
    `(progn
       (check-type ,numeric-vector numeric-vector)
       (check-type ,lla-type lla-type)
       (with-pinned-vector ((convert-elements ,numeric-vector ,lla-type))
         ,pointer)
       ,@body)))


#+sbcl
(defmacro with-nv-input-output ((numeric-vector output pointer lla-type)
				&body body)
  "A combination of with-nv-input and -output: input is always copied
\(and converted if necessary), is available during body, and the
contents end up in a _lisp vector_ of the appropriate type (ie it can
be used as elements in MAKE-NV* and MAKE-MATRIX*), assigned to name."
  (check-type output symbol)
  (check-type pointer symbol)
  (once-only (numeric-vector lla-type)
    `(progn
       (check-type ,numeric-vector numeric-vector)
       (check-type ,lla-type lla-type)
       (let* ((,output (convert-elements ,numeric-vector ,lla-type)))
	 (with-pinned-vector (,output ,pointer)
	   ,@body)))))

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
	 (with-pinned-vector ((elements ,name) ,pointer)
	   ,@body)))))
