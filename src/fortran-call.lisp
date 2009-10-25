(in-package :lla)

;;;; Code for interfacing with Fortran calls.

;;;; Ideally, this would not be the place for LAPACK-specific stuff,
;;;; just things one would generally use to call Fortran procedures.
;;;; However, for experimenting I will keep other idiosyncratic stuff
;;;; here.


;;;; ** LAPACK specific utility functions, macros and conditions.
;;;;
;;;;

(defmacro lb-procedure-name (name lla-type)
  "Expand to an (ecase lla-type) which should evaluate to the
procedure name.  LLA-TYPE has to be a symbol, if you need conditionals
etc, do that outside this macro."
  (check-type name symbol)
  (check-type lla-type symbol)
  `(ecase ,lla-type
     (:single ',(make-symbol* "%S" name))
     (:double ',(make-symbol* "%D" name))
     (:complex-single ',(make-symbol* "%C" name))
     (:complex-double ',(make-symbol* "%Z" name))))

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

(defmethod lu ((a dense-matrix-like))
  (as-dense-matrix (a)
    (bind (((:slots-read-only (m nrow) (n ncol) (a-data data)) a)
           (type (lla-type a-data))
           (procedure (lb-procedure-name getrf type)))
      (with-nv-input-output (a-data lu-data a% type)
        (with-nv-output (ipiv (min m n) ipiv% :integer)
          (with-fortran-scalars ((m m% :integer)
                                 (n n% :integer))
	    (with-info-check (getrf info%)
	      (funcall procedure m% n% a% m% ipiv% info%)))
          (make-instance 'lu :nrow m :ncol n :ipiv ipiv :data lu-data))))))

(defgeneric solve (a b)
  (:documentation "Return X that solves AX=B."))

(defmethod solve (a (b numeric-vector))
  ;; will simply call one of the methods below
  (matrix->vector (solve a (vector->matrix-col b))))

(defmethod solve ((a dense-matrix-like) (b dense-matrix-like))
  (as-dense-matrix (a b)
    (bind (((:slots-read-only (n nrow) (n2 ncol) (a-data data)) a)
           ((:slots-read-only (n3 nrow) (nrhs ncol) (b-data data)) b)
           (common-type (smallest-common-target-type 
                         (mapcar #'lla-type (list a-data b-data))))
           (procedure (lb-procedure-name gesv common-type)))
      (assert (= n n2 n3))
      (with-nv-input-copied (a-data a% common-type) ; LU decomposition discarded
        (with-nv-input-output (b-data x-data b% common-type)
          (with-work-area (ipiv% :integer n) ; not saving IPIV
            (with-fortran-scalars ((n n% :integer)
                                   (nrhs nrhs% :integer))
              (with-info-check (gesv info%)
                (funcall procedure n% nrhs% a% n% ipiv% b% n% info%)))
            (make-instance 'dense-matrix :nrow n :ncol nrhs
                           :data x-data)))))))

(defmethod solve ((lu lu) (b dense-matrix-like))
  (as-dense-matrix (b)
    (bind (((:slots-read-only (n nrow) (n2 ncol) ipiv (lu-data data)) lu)
           ((:slots-read-only (n3 nrow) (nrhs ncol) (b-data data)) b)
           (common-type (smallest-common-target-type
                         (mapcar #'lla-type (list lu-data b-data))))
           (procedure (lb-procedure-name getrs common-type)))
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
                             :data x-data))))))))

(defun eigen-dense-matrix-double (a &key vectors-p check-real-p)
  "Eigenvalues and vectors for dense, double matrices."
  (as-dense-matrix (a)
    (bind (((:slots-read-only (n nrow) (n2 ncol) (a-data data)) a)
           (procedure (lb-procedure-name geev :double)))
      (assert (= n n2))
      (with-nv-input-copied (a-data a% :double) ; A is overwritten
        (with-work-area (w% :double (* 2 n)) ; eigenvalues, will be zipped
          (let (;; imaginary part
                (wi% (inc-pointer w% (* n (foreign-type-size :double)))))
            (with-fortran-scalar (n n% :integer)
              (with-characters ((n-char #\N)
                                (v-char #\V))
                (if vectors-p
                    (with-nv-output (vr-data (* n n) vr% :double)
                      (with-lwork-query (lwork work :double)
                        (with-info-check (geev info%)
                          (funcall procedure n-char v-char n% a% n% 
                                   w% wi% ; eigenvalues
                                   (null-pointer) n% vr% n% ; eigenvectors
                                   work lwork info%)))
                      (values (nv-zip-complex-double w% n check-real-p)
                              (make-instance 'dense-matrix :nrow n :ncol n
                                             :data vr-data)))
                    (with-lwork-query (lwork work :double)
                      (with-info-check (geev info%)
                        (funcall procedure n-char n-char n% a% n% 
                                 w% wi% ; eigenvalues
                                 (null-pointer) n% (null-pointer) n% ; eigenvectors
                                 work lwork info%)
                        (nv-zip-complex-double w% n check-real-p))))))))))))



(defun eigen-dense-matrix-complex-double (a &key vectors-p check-real-p)
  "Eigenvalues and vectors for dense, double matrices."
  (declare (ignore a vectors-p check-real-p))
  (error "this function needs to be implemented"))


(defgeneric eigen (a &key vectors-p check-real-p &allow-other-keys)
  (:documentation "Calculate the eigenvalues and optionally the right
eigenvectors of a matrix.  Return (values eigenvalues eigenvectors).
If check-real-p, eigenvalues of real matrices are checked for an
imaginary part and returned with the appropriate type (compex or
not)."))

(defmethod eigen ((a dense-matrix-like) &key vectors-p check-real-p)
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

(defmethod least-squares ((a dense-matrix-like) (b dense-matrix-like))
  (as-dense-matrix (a b)
    (bind (((:slots-read-only (m nrow) (n ncol) (a-data data)) a)
           ((:slots-read-only (m2 nrow) (nrhs ncol) (b-data data)) b)
           (common-type (smallest-common-target-type
                         (mapcar #'lla-type (list a-data b-data))))
           (procedure (lb-procedure-name gels common-type)))
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
                  (funcall procedure n-char m% n% nrhs% a% m% b% m%
                           work% lwork% info%)))
              (values 
                (matrix-from-first-rows x-data m nrhs n)
                (make-instance 'qr :nrow m :ncol n
                               :data qr-data)
                (sum-last-rows x-data m nrhs n)))))))))

;;; univariate versions of least squares: vector ~ vector, vector ~ matrix

(defmethod least-squares ((a dense-matrix-like) (b numeric-vector))
  (bind (((:values x qr ss) (least-squares a (vector->matrix-col b))))
    (values (matrix->vector x) qr (xref ss 0))))

(defmethod least-squares ((a numeric-vector) (b numeric-vector))
  (bind (((:values x qr ss) (least-squares (vector->matrix-col a)
                                           (vector->matrix-col b))))
    (values (matrix->vector x) qr (xref ss 0))))


;;;;
;;;; inverting matrices
;;;;

(defgeneric invert (a)
  (:documentation "Invert A."))

(defmethod invert ((a dense-matrix-like))
  (as-dense-matrix (a)
    (invert (lu a))))

(defmethod invert ((lu lu))
  (bind (((:slots-read-only (n nrow) (n2 ncol) ipiv (lu-data data)) lu)
	 (common-type (lla-type lu-data))
	 (procedure (lb-procedure-name getri common-type)))
    (assert (= n n2))
    (with-nv-input-output (lu-data inverse-data lu% common-type)
      (with-nv-input (ipiv ipiv% :integer)
        (with-fortran-scalar (n n% :integer)
          (with-lwork-query (lwork% work% common-type)
            (with-info-check (getri info%)
		(funcall procedure n% lu% n% ipiv% work% lwork% info%)))
	  (make-instance 'dense-matrix :nrow n :ncol n
			 :data inverse-data))))))


(defun invert-triangular% (a upper-p unit-diag-p result-class)
  "Invert a dense (triangular) matrix using the LAPACK routine *TRTRI.
UPPER-P indicates if the matrix is in the upper or the lower triangle
of a (which needs to be a subtype of dense-matrix-like, but the type
information is not used), UNIT-DIAG-P indicates whether the diagonal
is supposed to consist of 1s.

NOTE: for internal use, not exported."
  (bind (((:slots-read-only (n nrow) (n2 ncol) (a-data data)) a)
	 (common-type (lla-type a-data))
	 (procedure (lb-procedure-name trtri common-type)))
    (assert (= n n2))
    (with-nv-input-output (a-data inv-data a% common-type)
      (with-characters ((u-char (if upper-p #\U #\L))
                        (d-char (if unit-diag-p #\U #\N)))
        (with-fortran-scalar (n n% :integer)
          (with-info-check (trtri info%)
            (funcall procedure u-char d-char n% a% n% info%))))
      (make-instance result-class :nrow n :ncol n :data inv-data))))
  
(defmethod invert ((a upper-triangular-matrix))
  (invert-triangular% a t nil 'upper-triangular-matrix))

(defmethod invert ((a lower-triangular-matrix))
  (invert-triangular% a nil nil 'lower-triangular-matrix))

(defmethod invert ((a cholesky))
  (bind (((:slots-read-only (n nrow) (n2 ncol) (a-data data)) a)
	 (common-type (lla-type a-data))
	 (procedure (lb-procedure-name potri common-type)))
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
;;;;  ??? In the future, C and beta could be added.  Is it needed? --
;;;;  Tamas

(defgeneric mm (a b &key alpha)
  (:documentation "multiply A and B, also by the scalar
  alpha (defaults to 1)."))

(defmethod mm ((a dense-matrix-like) (b dense-matrix-like) &key (alpha 1))
  (as-dense-matrix (a b)
    (bind (((:slots-read-only (m nrow) (k ncol) (a-data data)) a)
           ((:slots-read-only (k2 nrow) (n ncol) (b-data data)) b)
           (common-type (smallest-common-target-type
                         (mapcar #'lla-type (list a-data b-data))))
           (lisp-type (lla-type->lisp-type common-type))
           (procedure (lb-procedure-name gemm common-type)))
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
            (make-instance 'dense-matrix :nrow m :ncol n :data c-data)))))))

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
;; 	 (procedure (lb-procedure-name 'trmm common-type)))
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

(defmethod mm ((a numeric-vector) (b dense-matrix-like) &key (alpha 1))
  (matrix->vector (mm (vector->matrix-row a) b :alpha alpha)))

(defmethod mm ((a dense-matrix-like) (b numeric-vector) &key (alpha 1))
  (matrix->vector (mm a (vector->matrix-col b) :alpha alpha)))

(defmethod mm ((a numeric-vector) (b numeric-vector) &key (alpha 1))
  ;; a dot product, basically :-)
  (xref (mm (vector->matrix-row a) (vector->matrix-col b) :alpha alpha)
        0 0))
  




;;;;
;;;;  Updates for symmetric and hermitian matrices
;;;;

(defun ntc-operation-char% (operation)
  "Return corresponding character for LAPACK operations."
  (ecase operation
    (:none #\N)
    (:transpose #\T)
    (:conjugate #\C)))

(defgeneric update-syhe (a c op-left-p &key alpha beta)
  (:documentation "Calculate alpha*A*op(A) + beta*C if OP-LEFT-P is
  NIL, and alpha*op(A)*A + beta*C otherwise.  ALPHA and BETA default
  to 1 and 0, respectively.

  If C is a HERMITIAN-MATRIX or equal to 'HERMITIAN-MATRIX (in which
  case C is created as such, filled with 0s), then op(A) is conjugate
  transpose, and the result is a HERMITIAN-MATRIX.

  If C is a SYMMETRIC-MATRIX or equal to 'SYMMETRIC-MATRIX (in which
  case C is created as such, filled with 0s), then op(A) is transpose,
  and the result is a SYMMETRIC-MATRIX."))

(defun update-syhe% (a c alpha beta operation hermitian-p)
  "Internal procedure for *SYRK and *HERK.  The convention is that
this function will convert A to a DENSE-MATRIX."
  (as-dense-matrix (a)
    (bind (((:slots-read-only (nrow-a nrow) (ncol-a ncol) (a-data data)) a)
           ((:slots-read-only (nrow-c nrow) (ncol-c ncol) (c-data data)) c)
           (common-type (let ((ct (smallest-common-target-type
                                   (mapcar #'lla-type (list c-data c-data)))))
                          (if hermitian-p
                              ;; for Hermitian operations, the result is
                              ;; always complex.
                              (case ct
                                ((:single :complex-single) :complex-single)
                                (t :complex-double))
                              ct)))
           (lisp-type (lla-type->lisp-type common-type))
           (procedure (if hermitian-p
                          (lb-procedure-name herk common-type)
                          (lb-procedure-name syrk common-type)))
           (k (ecase operation          ; the "other" dimension of A
                (:none (assert (= nrow-a nrow-c)) ncol-a)
                ((:transpose :conjugate) (assert (= ncol-a ncol-c)) nrow-a))))
      (assert (= ncol-c nrow-c)) ; not really needed, c has to be symm/herm
      (when hermitian-p
        ;; T not a valid operation, LAPACK would choke
        (assert (not (eq operation :transpose))))
      (with-nv-input (a-data a% common-type)
        (with-nv-input-output (c-data x-data c% common-type)
          (with-fortran-scalars ((nrow-a nrow-a% :integer)
                                 (ncol-c n% :integer)
                                 (k k% :integer)
                                 ((coerce alpha lisp-type) alpha% common-type)
                                 ((coerce beta lisp-type) beta% common-type))
            (with-characters ((u-char #\U)
                              (t-char (ntc-operation-char% operation)))
              (funcall procedure u-char t-char n% k% alpha% a% nrow-a%
                       beta% c% n%)))
          (make-instance (if hermitian-p
                             'hermitian-matrix
                             'symmetric-matrix) 
                         :nrow nrow-c :ncol nrow-c :data x-data))))))

(defmethod update-syhe ((a dense-matrix-like) (c symmetric-matrix) op-left-p
                        &key (alpha 1) (beta 0))
    (update-syhe% a c alpha beta (if op-left-p :transpose :none) nil))

(defmethod update-syhe ((a dense-matrix-like) (c hermitian-matrix) op-left-p
                        &key (alpha 1) (beta 0))
  (let ((operation (if op-left-p :conjugate :none)))
    (if (or (lla-complex-p (lla-type a)) (lla-complex-p (lla-type c)))
        (update-syhe% a c alpha beta operation t)
        (take 'hermitian-matrix
              (update-syhe% a c alpha beta operation nil)))))

(defun update-syhe-creating-c% (a c op-left-p alpha beta)
  "Create a conforming C matrix and clal update-syhe.  C has to be one
of the valid symbols.  Internal, not exported."
  (unless (zerop beta)
    (warn "beta=~a does not make much sense for a zero c" beta))
  (let* ((n (xdim a (if op-left-p 1 0)))
         (c (make-matrix c n n :lla-type (lla-type a) :initial-contents 0)))
    (update-syhe a c op-left-p :alpha alpha :beta beta)))

(defmethod update-syhe ((a dense-matrix-like) (c (eql 'symmetric-matrix)) op-left-p
                        &key (alpha 1) (beta 0))
  (update-syhe-creating-c% a c op-left-p alpha beta))

(defmethod update-syhe ((a dense-matrix-like) (c (eql 'hermitian-matrix)) op-left-p
                        &key (alpha 1) (beta 0))
  (update-syhe-creating-c% a c op-left-p alpha beta))


;;;;
;;;;  Cholesky factorization
;;;;

(defgeneric cholesky (a)
  (:documentation "Cholesky factorization.  Only uses the lower
  triangle of a dense-matrix, and needs a PSD matrix."))

(defmethod cholesky ((a dense-matrix-like))
  (as-dense-matrix (a)
    (bind (((:slots-read-only (n nrow) (n2 ncol) data) a)
           (common-type (lla-type data))
           (procedure (lb-procedure-name potrf common-type)))
      (assert (= n n2))
      (with-nv-input-output (data cholesky-data data% common-type)
        (with-fortran-scalar (n n% :integer)
          (with-character (u-char #\L)
            (with-info-check (potrf info%)
              (funcall procedure u-char n% data% n% info%)))
	  (make-instance 'cholesky :nrow n :ncol n
			 :data cholesky-data))))))

(defmethod reconstruct ((mf cholesky))
  (update-syhe (take 'lower-triangular-matrix mf) 'hermitian-matrix nil))
