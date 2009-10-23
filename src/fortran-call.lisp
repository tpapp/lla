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
             ;; ??? two different conditions should be thrown,
             ;; depending on the value of INFO.  Positive INFO usually
             ;; means something substantive (eg a minor is not PSD),
             ;; negative INFO is a bad argument which should never
             ;; happen
             (error 'lapack-error :info ,info-value 
                    :lapack-procedure ',procedure-name)))))))

(defmacro with-lwork-query ((lwork work lla-type) &body body)
  "Call body twice with the given parameters, querying the size for
the workspace area.  NOTE: abstraction leaks a bit (body is there
twice), but it should not be a problem in practice.  Body is most
commonly a single function call."
  (check-type lwork symbol)
  (check-type work symbol)
  (with-unique-names (work-size foreign-size)
    (once-only (lla-type)
      `(with-foreign-object (,lwork :int32 1)
         (let ((,foreign-size (foreign-size* ,lla-type))
               ,work-size)
           ;; query workspace size
           (setf (mem-ref ,lwork :int32) -1)
           (with-foreign-pointer (,work ,foreign-size)
             ,@body
             (setf ,work-size (floor (mem-aref* ,work ,lla-type))
                   (mem-ref ,lwork :int32) ,work-size))
           ;; allocate and use workspace
           (with-foreign-pointer (,work (* ,foreign-size ,work-size))
             ,@body))))))


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
	 (type (lla-type a-data))
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
	 (common-type (smallest-common-target-type 
                       (mapcar #'lla-type (list a-data b-data))))
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
	 (common-type (smallest-common-target-type
                       (mapcar #'lla-type (list lu-data b-data))))
	 (procedure (lapack-procedure-name 'getrs common-type)))
    (assert (= n n2 n3))
    (with-nv-input (lu-data lu% common-type)
      (with-nv-input-output (b-data x-data b% common-type)
	(with-nv-input (ipiv ipiv% :integer)
	  (with-fortran-scalars ((n n% :integer)
				 (nrhs nrhs% :integer))
	    (with-character (trans% #\N)
	      (with-info-check (getrf info%)
		(funcall procedure trans% n% nrhs% lu% n% ipiv% b% n% info%)))
	  (make-instance 'dense-matrix :nrow n :ncol nrhs
			 :data x-data)))))))

(defun eigen-dense-matrix-double (a &key vectors-p check-real-p)
  "Eigenvalues and vectors for dense, double matrices."
  (bind (((:slots-read-only (n nrow) (n2 ncol) (a-data data)) a)
	 (procedure (lapack-procedure-name 'geev :double)))
    (assert (= n n2))
    (with-nv-input-copied (a-data a% :double) ; A is overwritten
      (with-work-area (w% :double (* 2 n))     ; eigenvalues, will be zipped
        (let ((wi% (inc-pointer w% (* n (foreign-type-size :double))))) ; imaginary part
          (with-fortran-scalar (n n% :integer)
            (with-characters ((n-char #\N)
                              (v-char #\V))
              (if vectors-p
                  (with-nv-output (vr-data (* n n) vr% :double)
                    (with-lwork-query (lwork work :double)
                      (with-info-check (geev info%)
                        (funcall procedure n-char v-char n% a% n% 
                                 w% wi%                   ; eigenvalues
                                 (null-pointer) n% vr% n% ; eigenvectors
                                 work lwork info%)))
                    (values (nv-zip-complex-double w% n check-real-p)
                            (make-instance 'dense-matrix :nrow n :ncol n
                                           :data vr-data)))
                  (with-lwork-query (lwork work :double)
                    (with-info-check (geev info%)
                      (funcall procedure n-char n-char n% a% n% 
                             w% wi%                              ; eigenvalues
                             (null-pointer) n% (null-pointer) n% ; eigenvectors
                             work lwork info%)
                      (nv-zip-complex-double w% n check-real-p)))))))))))



(defun eigen-dense-matrix-complex-double (a &key vectors-p check-real-p)
  "Eigenvalues and vectors for dense, double matrices."
  (declare (ignore a vectors-p check-real-p))
  (error "this function needs to be written"))


(defgeneric eigen (a &key vectors-p check-real-p &allow-other-keys)
  (:documentation "Calculate the eigenvalues and optionally the right
eigenvectors of a matrix.  Return (values eigenvalues eigenvectors).
If check-real-p, eigenvalues of real matrices are checked for an
imaginary part and returned with the appropriate type (compex or
not)."))

(defmethod eigen ((a dense-matrix) &key vectors-p check-real-p)
  ;; The current approach is: convert to double precision (complex or
  ;; real), so we just need two functions.  Unfortunately, the LAPACK
  ;; interface for real and complex cases is different.
  (case (lla-type (data a))
    ((:integer :single :double)
       (eigen-dense-matrix-double a :vectors-p vectors-p
                                  :check-real-p check-real-p))
    ((:complex-single :complex-double)
       (eigen-dense-matrix-complex-double a :vectors-p vectors-p
                                          :check-real-p check-real-p))))

;;;;
;;;; least squares calculations
;;;;

(defgeneric least-squares (a b)
  (:documentation "Return (values x qr ss), where x = argmin_x L2norm(
b-Ax ), solving a least squares problem, qr is the QR decomposition of
A, and SS is the sum of squares for each column of B.  B can have
multiple columns, in which case x will have the same number of
columns, each corresponding to a different column of b."))

(defmethod least-squares ((a dense-matrix) (b dense-matrix))
  (bind (((:slots-read-only (m nrow) (n ncol) (a-data data)) a)
	 ((:slots-read-only (m2 nrow) (nrhs ncol) (b-data data)) b)
	 (common-type (smallest-common-target-type
                       (mapcar #'lla-type (list a-data b-data))))
	 (procedure (lapack-procedure-name 'gels common-type)))
    (assert (= m m2))
    (unless (<= n m)
      (error "A doesn't have enough columns for least squares"))
    (with-nv-input-output (a-data qr-data a% common-type) ; output: QR decomposition
      (with-nv-input-output (b-data x-data b% common-type) ; output: 
        (with-character (n-char #\N)
	  (with-fortran-scalars ((n n% :integer)
                                 (m m% :integer)
				 (nrhs nrhs% :integer))
            (with-lwork-query (lwork% work% :double)
              (with-info-check (gels info%)
                (funcall procedure n-char m% n% nrhs% a% m% b% m% work% lwork% info%)))
            (values 
              (matrix-from-first-rows x-data m nrhs n)
              (make-instance 'qr :nrow m :ncol n
                             :data qr-data)
              (sum-last-rows x-data m nrhs n))))))))

;;; univariate versions of least squares: vector ~ vector, vector ~ matrix

(defmethod least-squares ((a dense-matrix) (b numeric-vector))
  (bind (((:values x qr ss) (least-squares a (vector->matrix-col b))))
    (values (matrix->vector x) qr (xref ss 0))))

(defmethod least-squares ((a numeric-vector) (b numeric-vector))
  (bind (((:values x qr ss) (least-squares (vector->matrix-col a) (vector->matrix-col b))))
    (values (matrix->vector x) qr (xref ss 0))))


;;;;
;;;; inverting matrices
;;;;

(defgeneric invert (a)
  (:documentation "Invert A."))

(defmethod invert ((a dense-matrix))
  (invert (lu a)))

(defmethod invert ((lu lu))
  (bind (((:slots-read-only (n nrow) (n2 ncol) ipiv (lu-data data)) lu)
	 (common-type (lla-type lu-data))
	 (procedure (lapack-procedure-name 'getri common-type)))
    (assert (= n n2))
    (with-nv-input-output (lu-data inverse-data lu% common-type)
      (with-nv-input (ipiv ipiv% :integer)
        (with-fortran-scalar (n n% :integer)
          (with-lwork-query (lwork% work% common-type)
            (with-info-check (getri info%)
		(funcall procedure n% lu% n% ipiv% work% lwork% info%)))
	  (make-instance 'dense-matrix :nrow n :ncol n
			 :data inverse-data))))))


(defun invert-triangular (a upper-p unit-diag-p result-class)
  "Invert a dense (triangular) matrix using the LAPACK routine *TRTRI.
UPPER-P indicates if the matrix is in the upper or the lower triangle
of a (which needs to be a subtype of dense-matrix, but the type
information is not used), UNIT-DIAG-P indicates whether the diagonal
is supposed to consist of 1s.

NOTE: for internal use, not exported."
  (bind (((:slots-read-only (n nrow) (n2 ncol) (a-data data)) a)
	 (common-type (lla-type a-data))
	 (procedure (lapack-procedure-name 'trtri common-type)))
    (assert (= n n2))
    (with-nv-input-output (a-data inv-data a% common-type)
      (with-characters ((u-char (if upper-p #\U #\L))
                        (d-char (if unit-diag-p #\U #\N)))
        (with-fortran-scalar (n n% :integer)
          (with-info-check (trtri info%)
            (funcall procedure u-char d-char n% a% n% info%))))
      (make-instance result-class :nrow n :ncol n :data inv-data))))
  
(defmethod invert ((a upper-triangular-matrix))
  (invert-triangular a t nil 'upper-triangular-matrix))

(defmethod invert ((a lower-triangular-matrix))
  (invert-triangular a nil nil 'lower-triangular-matrix))

(defmethod invert ((a cholesky))
  (bind (((:slots-read-only (n nrow) (n2 ncol) (a-data data)) a)
	 (common-type (lla-type a-data))
	 (procedure (lapack-procedure-name 'potri common-type)))
    (assert (= n n2))
    (with-nv-input-output (a-data inv-data a% common-type) ; output: upper triangular
      (with-character (u-char #\U)
        (with-fortran-scalar (n n% :integer)
          (with-info-check (gels info%)
            (funcall procedure u-char n% a% n% info%))))
      (make-instance 'symmetric-matrix :nrow n :ncol n :data inv-data)))) 

;;;;
;;;;  utility functions for least squares
;;;;

(defun least-squares-raw-variance (qr)
  "Calculate the residual variance (basically, (X^T X)-1 ) from the qr
decomposition of X."
  ;; Notes: X = QR, thus X^T X = R^T Q^T Q R = R^T R because Q is
  ;; orthogonal.  Then we do as if calculating the inverse of a matrix
  ;; using its Cholesky decomposition.
  (with-slots (nrow ncol data) qr
    (assert (<= ncol nrow))
    (invert (take 'cholesky (factorization-component qr :R)))))


;;;;
;;;;  matrix multiplication
;;;; 
;;;;  Currently, we only support the common alpha parameter alpha*A*B.

(defgeneric mm (a b &optional alpha)
  (:documentation "multiply A and B"))

(defmethod mm ((a dense-matrix) (b dense-matrix) &optional (alpha 1))
  (bind ((a (take 'dense-matrix a))
         (b (take 'dense-matrix b))
         ((:slots-read-only (m nrow) (k ncol) (a-data data)) a)
	 ((:slots-read-only (k2 nrow) (n ncol) (b-data data)) b)
	 (common-type (smallest-common-target-type
		       (mapcar #'lla-type (list a-data b-data))))
         (lisp-type (lla-type->lisp-type common-type))
	 (procedure (lapack-procedure-name 'gemm common-type)))
    (assert (= k k2))
    (with-nv-input (a-data a% common-type)
      (with-nv-input (b-data b% common-type)
	(with-nv-output (c-data (* m n) c% common-type)
	  (with-fortran-scalars ((n n% :integer)
				 (k k% :integer)
				 (m m% :integer)
				 ((coerce alpha lisp-type) alpha% common-type)
				 ((coerce 0 lisp-type) z% common-type))
	    (with-character (trans% #\N)
	      (funcall procedure trans% trans% m% n% k% alpha% a% m% b% k%
		       z% c% m%)))
	  (make-instance 'dense-matrix :nrow m :ncol n :data c-data))))))

;;; !!! write triangular method, currently I am confused about its
;;; !!! generality wrt dimensions - does it handle square triangular
;;; !!! matrices only? currently the dense-matrix version will do just fine.

;; (defun mm-triangular (a b side upper-p unit-diag-p alpha)
;;   "Utility function for multiplication of triangular matrices.  SIDE
;; is either :LEFT (AB) or :RIGHT (BA), but A is always the triangular
;; matrix.  NOTE: not exported, for internal use only."
;;   (bind (((:slots-read-only (nrow-a nrow) (ncol-a ncol) (a-data data)) a)
;; 	 ((:slots-read-only (nrow-b nrow) (ncol-b ncol) (b-data data)) b)
;; 	 (common-type (smallest-common-target-type
;; 		       (mapcar #'lla-type (list a-data b-data))))
;;          (lisp-type (lla-type->lisp-type common-type))
;; 	 (procedure (lapack-procedure-name 'trmm common-type)))
;;     (case side
;;       (:left (assert (= ncol-a nrow-b)))
;;       (:right (assert (= ncol-b nrow-a)))
;;       (otherwise (error "invalide side specification ~A" side)))
;;     (with-nv-input (a-data a% common-type)
;;       (with-nv-input-output (b-data x-data b% common-type)
;;         (with-fortran-scalars ((nrow-a nrow-a% :integer)
;;                                (ncol-a ncol-a% :integer)
;;                                (nrow-b nrow-b% :integer)
;;                                (ncol-b ncol-b% :integer)
;;                                ((coerce alpha lisp-type) alpha% common-type)
;;                                ((coerce 0 lisp-type) z% common-type))
;;           (with-characters ((t-char #\N)
;;                             (u-char (if upper-p #\U #\L))
;;                             (d-char (if unit-diag-p #\U #\N))
;;                             (s-char (if (eq side :left) #\L #\R)))
;;             (funcall procedure s-char u-char t-char d-char
;;                      nrow-b% ncol-b% alpha% a% )))
;; 	  (make-instance 'dense-matrix :nrow m :ncol n :data c-data))))))


;; (defmethod mm ((a upper-triangular-

(defmethod mm ((a numeric-vector) (b dense-matrix) &optional (alpha 1))
  (matrix->vector (mm (vector->matrix-row a) b alpha)))

(defmethod mm ((a dense-matrix) (b numeric-vector) &optional (alpha 1))
  (matrix->vector (mm a (vector->matrix-col b) alpha)))

;;;;
;;;;  Cholesky factorization
;;;;

(defgeneric cholesky (a)
  (:documentation "Cholesky factorization.  Only used the upper
  triangle of a dense-matrix, and needs a PSD matrix."))

(defmethod cholesky ((a dense-matrix))
  (bind (((:slots-read-only (n nrow) (n2 ncol) data) a)
	 (common-type (lla-type data))
	 (procedure (lapack-procedure-name 'potrf common-type)))
    (assert (= n n2))
    (with-nv-input-output (data cholesky-data data% common-type)
      (with-fortran-scalar (n n% :integer)
        (with-character (u-char #\U)
          (with-info-check (potrf info%)
            (funcall procedure u-char n% data% n% info%)))
	  (make-instance 'cholesky :nrow n :ncol n
			 :data cholesky-data)))))

(defmethod reconstruct ((mf cholesky))
  (let ((herm (take 'hermitian-matrix mf)))
    (set-restricted herm)              ; mirror on the diagonal
    (let* ((ut (take 'upper-triangular-matrix herm))
           (lt (take 'lower-triangular-matrix herm))
           (result (mm lt ut)))
      (take (if (lla-complex-p (lla-type result))
                'hermitian-matrix
                'symmetric-matrix)
            result))))
