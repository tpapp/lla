(in-package :lla)

;;;; Code for interfacing with Fortran calls.

;;;; Ideally, this would not be the place for LAPACK-specific stuff,
;;;; just things one would generally use to call Fortran procedures.
;;;; However, for experimenting I will keep other idiosyncratic stuff
;;;; here.


;;;; ** LAPACK specific utility functions, macros and conditions.
;;;;
;;;;

(defun lapack-procedure-name (name lla-type)
  "Return the name for the LAPACK procedure name based on name and
general-type.  Name can be a string or a symbol, the returned value is
a symbol."
  (make-symbol* (ecase lla-type
		  (:single "%S")
		  (:double "%D")
		  (:complex-single "%C")
		  (:complex-double "%Z"))
		name))

(define-condition lapack-error (error)
  ;; !! write method for formatting the error message
  ((info :initarg :info :type integer :reader info
	 :documentation "info code")
   (lapack-procedure :initarg :lapack-procedure :type symbol
		     :reader lapack-procedure
		     :documentation "The name without the type prefix (eg 'gesv)."))
  (:documentation "The lapack procedure returned a nonzero info
  code."))

(defmacro with-info-check ((procedure-name info-pointer) &body body)
  "Evaluate body with info-pointer bound to an integer, and check it
afterwards, signalling an lapack-error condition if info is nonzero."
  ;;;; !!! how to handle errors nicely? the wrapper can handle this
  ;;;; condition, and output an error message.  There are generally
  ;;;; two kinds of errors in LAPACK: (1) malformed inputs, and (2)
  ;;;; ill-conditioned numerical problems (eg trying to invert a
  ;;;; matrix which is not invertible, etc).
  ;;;;
  ;;;; (1) requires displaying argument names, but in a well-written
  ;;;; library errors like that should not happen, so we will not
  ;;;; worry about it (if it does, that requires debugging & fixing,
  ;;;; not a condition system). In (2), info usually carries all the
  ;;;; information.
  ;;;;
  ;;;; ??? suggestion: an extra argument to the macro on how to
  ;;;; interpret info and generate an error message? create a class
  ;;;; hierarchy of conditions. !!!! do this when API has
  ;;;; stabilized. -- Tamas
  (check-type info-pointer symbol)
  (check-type procedure-name symbol)
  (with-unique-names (info-value)
    `(with-foreign-object (,info-pointer :int32)
       (multiple-value-prog1
	   (progn ,@body)
	 (let ((,info-value (mem-aref ,info-pointer :int32)))
	   (unless (zerop ,info-value)
		 (error 'lapack-error :info ,info-value 
			:lapack-procedure ',procedure-name)))))))

;;;; Actual functions defined here.  Eventually, they will be moved to
;;;; somewhere else.  And more importantly, I am planning to write a
;;;; macro system to generate these functions intelligently from a few
;;;; lines of specs (see notused/utilities.lisp for previous version),
;;;; as there is a lot of repetitive stuff going on.  This is just to
;;;; explore what that macro needs to generate.

;;;; General convention for vectors in places of matrices: should be
;;;; interpreted as a conforming vector.  If the result is an 1xn or
;;;; nx1 matrix, it should be converted to a vector iff some related
;;;; argument was a vector.  For example, (solve a b) should be a
;;;; vector iff b is a vector, otherwise it should be a matrix.

(defgeneric lu (a)
  (:documentation "LU decomposition of A"))

(defmethod lu ((a dense-matrix))
  (bind (((:slots-read-only (m nrow) (n ncol) (a-data data)) a)
	 (type (element-lla-type a-data))
	 (procedure (lapack-procedure-name 'getrf type)))
    (with-nv-input-output (a-data lu-data a% type)
      (with-nv-output (ipiv (min m n) ipiv% :integer)
	(with-fortran-scalars ((m m% :integer)
			       (n n% :integer))
	    (with-info-check (getrf info%)
	      (funcall procedure m% n% a% m% ipiv% info%)))
	(make-instance 'lu :nrow m :ncol n :ipiv ipiv :data lu-data)))))

(defgeneric solve (a b)
  (:documentation "Return X that solves AX=B."))

(defmethod solve (a (b numeric-vector))
  ;; will simply call one of the methods below
  (matrix->vector (solve a (vector->matrix-col b))))

(defmethod solve ((a dense-matrix) (b dense-matrix))
  (bind (((:slots-read-only (n nrow) (n2 ncol) (a-data data)) a)
	 ((:slots-read-only (n3 nrow) (nrhs ncol) (b-data data)) b)
	 (common-type (smallest-common-target-type (mapcar #'element-lla-type (list a-data b-data))))
	 (procedure (lapack-procedure-name 'gesv common-type)))
    (assert (= n n2 n3))
    (with-nv-input-copied (a-data a% common-type) ; LU decomposition discarded
      (with-nv-input-output (b-data x-data b% common-type)
	(with-work-area (ipiv% :integer n) ; not saving IPIV
	  (with-fortran-scalars ((n n% :integer)
				 (nrhs nrhs% :integer))
	    (with-info-check (gesv info%)
	      (funcall procedure n% nrhs% a% n% ipiv% b% n% info%)))
	  (make-instance 'dense-matrix :nrow n :ncol nrhs
			 :data x-data))))))

(defmethod solve ((lu lu) (b dense-matrix))
  (bind (((:slots-read-only (n nrow) (n2 ncol) ipiv (lu-data data)) lu)
	 ((:slots-read-only (n3 nrow) (nrhs ncol) (b-data data)) b)
	 (common-type (smallest-common-target-type (mapcar #'element-lla-type (list lu-data b-data))))
	 (procedure (lapack-procedure-name 'getrs common-type)))
    (assert (= n n2 n3))
    (with-nv-input (lu-data lu% common-type)
      (with-nv-input-output (b-data x-data b% common-type)
	(with-nv-input (ipiv ipiv% :integer)
	  (with-fortran-scalars ((n n% :integer)
				 (nrhs nrhs% :integer))
	    (with-characters ((#\N trans%))
	      (with-info-check (getrf info%)
		(funcall procedure trans% n% nrhs% lu% n% ipiv% b% n% info%)))
	  (make-instance 'dense-matrix :nrow n :ncol nrhs
			 :data x-data)))))))

(defgeneric mm (a b)
  (:documentation "multiply A and B"))

(defmethod mm ((a dense-matrix) (b dense-matrix))
  (bind (((:slots-read-only (m nrow) (k ncol) (a-data data)) a)
	 ((:slots-read-only (k2 nrow) (n ncol) (b-data data)) b)
	 (common-type (smallest-common-target-type
		       (mapcar #'element-lla-type (list a-data b-data))))
	 (procedure (lapack-procedure-name 'gemm common-type)))
    (assert (= k k2))
    (with-nv-input (a-data a% common-type)
      (with-nv-input (b-data b% common-type)
	(with-nv-output (c-data (* m n) c% common-type)
	  (with-fortran-scalars ((n n% :integer)
				 (k k% :integer)
				 (m m% :integer)
				 (1d0 u% :double)
				 (0d0 z% :double))
	    (with-characters ((#\N trans%))
	      (funcall procedure trans% trans% m% n% k% u% a% m% b% k%
		       z% c% m%)))
	  (make-instance 'dense-matrix :nrow m :ncol n :data c-data))))))
