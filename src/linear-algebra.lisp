(in-package :lla)

;;;; Higher level linear algebra functions defined here.
;;;;
;;;; Eventually, I am planning to write a macro system to generate
;;;; these functions intelligently from a few lines of specs (see
;;;; notused/utilities.lisp for previous version), as there is a lot
;;;; of repetitive stuff going on.  This is just to explore what that
;;;; macro needs to generate.

;;;; General convention for vectors in places of matrices: should be
;;;; interpreted as a conforming vector.  If the result is an 1xn or
;;;; nx1 matrix, it should be converted to a vector iff some related
;;;; argument was a vector.  For example, (solve a b) should be a
;;;; vector iff b is a vector, otherwise it should be a matrix.

(defgeneric lu (a)
  (:documentation "LU decomposition of A"))

(defmethod lu ((a dense-matrix-like))
  (set-restricted a)
  (bind (((:slots-read-only (m nrow) (n ncol) (a-data data)) a)
         (type (lla-type a-data))
         (procedure (lb-procedure-name 'getrf type)))
    (with-nv-input-output (a-data lu-data a% type)
      (with-nv-output (ipiv (min m n) ipiv% :integer)
        (with-fortran-atoms ((:integer m% m)
                             (:integer n% n))
          (call-with-info-check procedure m% n% a% m% ipiv% info%))
        (make-instance 'lu :nrow m :ncol n :ipiv ipiv :data lu-data)))))

(defgeneric solve (a b)
  (:documentation "Return X that solves AX=B."))

(defmethod solve (a (b numeric-vector))
  ;; will simply call one of the methods below
  (matrix->vector (solve a (vector->matrix-col b))))

(defmethod solve ((a dense-matrix-like) (b dense-matrix-like))
  (set-restricted a)
  (set-restricted b)
  (bind (((:slots-read-only (n nrow) (n2 ncol) (a-data data)) a)
         ((:slots-read-only (n3 nrow) (nrhs ncol) (b-data data)) b)
         (common-type (smallest-common-target-type 
                       (mapcar #'lla-type (list a-data b-data))))
         (procedure (lb-procedure-name 'gesv common-type)))
    (assert (= n n2 n3))
    (with-nv-input-copied (a-data a% common-type) ; LU decomposition discarded
      (with-nv-input-output (b-data x-data b% common-type)
        (with-work-area (ipiv% :integer n) ; not saving IPIV
          (with-fortran-atoms ((:integer n% n)
                               (:integer nrhs% nrhs))
            (call-with-info-check procedure n% nrhs% a% n% ipiv% b% n% info%))
          (make-instance 'dense-matrix :nrow n :ncol nrhs
                         :data x-data))))))

(defmethod solve ((lu lu) (b dense-matrix-like))
  (set-restricted b)
  (bind (((:slots-read-only (n nrow) (n2 ncol) ipiv (lu-data data)) lu)
         ((:slots-read-only (n3 nrow) (nrhs ncol) (b-data data)) b)
         (common-type (smallest-common-target-type
                       (mapcar #'lla-type (list lu-data b-data))))
         (procedure (lb-procedure-name 'getrs common-type)))
    (assert (= n n2 n3))
    (with-nv-input (lu-data lu% common-type)
      (with-nv-input-output (b-data x-data b% common-type)
        (with-nv-input (ipiv ipiv% :integer)
          (with-fortran-atoms ((:integer n% n)
                               (:integer nrhs% nrhs)
                               (:char trans% #\N))
              (call-with-info-check procedure trans% n% nrhs% lu% n% ipiv% b% n% info%)))
        (make-instance 'dense-matrix :nrow n :ncol nrhs
                       :data x-data)))))

(defun eigen-dense-matrix-double (a &key vectors-p check-real-p)
  "Eigenvalues and vectors for dense, double matrices."
  (set-restricted a)
  (bind (((:slots-read-only (n nrow) (n2 ncol) (a-data data)) a)
         (procedure (lb-procedure-name 'geev :double)))
    (assert (= n n2))
    (with-nv-input-copied (a-data a% :double) ; A is overwritten
      (with-work-area (w% :double (* 2 n)) ; eigenvalues, will be zipped
        (let (;; imaginary part
              (wi% (inc-pointer w% (* n (foreign-type-size :double)))))
          (with-fortran-atoms ((:integer n% n)
                               (:char n-char #\N)
                               (:char v-char #\V))
              (if vectors-p
                  (with-nv-output (vr-data (* n n) vr% :double)
                    (with-work-query (lwork work :double)
                      (call-with-info-check procedure n-char v-char n% a% n% 
                                            w% wi% ; eigenvalues
                                            (null-pointer) n% vr% n% ; eigenvectors
                                            work lwork info%))
                    (values (nv-zip-complex-double w% n check-real-p)
                            (make-instance 'dense-matrix :nrow n :ncol n
                                           :data vr-data)))
                  (with-work-query (lwork work :double)
                    (call-with-info-check procedure n-char n-char n% a% n% 
                                          w% wi%   ; eigenvalues
                                          (null-pointer) n% (null-pointer) n% ; eigenvectors
                                          work lwork info%)
                    (nv-zip-complex-double w% n check-real-p)))))))))

(defun eigen-dense-matrix-complex-double (a &key vectors-p check-real-p)
  "Eigenvalues and vectors for dense, double matrices."
  (declare (ignore a vectors-p check-real-p))
  (error "this function needs to be implemented"))

(defgeneric eigen (a &key vectors-p check-real-p &allow-other-keys)
  (:documentation "Calculate the eigenvalues and optionally the right
eigenvectors of a matrix (as columns).  Return (values eigenvalues
eigenvectors).  If check-real-p, eigenvalues of real matrices are
checked for an imaginary part and returned with the appropriate
type (compex or not).  Complex conjugate pairs of eigenvalues appear
consecutively with the eigenvalue having the positive imaginary part
first."))

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

;;;  Currently, the setup below triggers floating point exceptions
;;;  (divison-by-zero), which is actually expected since it calls
;;;  *STEMR.  I plan to fix this with implementation-dependent flags
;;;  later on, currently we just use the simplified version.
;;;
;; (defmethod eigen ((a hermitian-matrix) &key vectors-p check-real-p)
;;   ;;;; Currently, we use the RRR (relatively robust representation)
;;   ;;;; algorithm for symmetric/hermitian eigenproblems.
;;   (declare (ignore check-real-p))       ; eigenvalues are always real
;;   (bind (((:slots-read-only (n nrow) (n2 ncol) (a-data data)) a)
;;          (lla-type (lla-type a-data))
;;          ((:values procedure real-lla-type) (lb-procedure-name2 'syevr 'heevr lla-type))
;;          (null% (null-pointer)))
;;     (assert (= n n2))
;;     (with-nv-input-copied (a-data a% lla-type) ; A is overwritten
;;       (with-nv-output (w n w% real-lla-type)   ; eigenvalues
;;         (with-fortran-atoms ((:integer n% n)
;;                              (real-lla-type abstol% (coerce* 0 real-lla-type)) ; function will select
;;                              (:char n-char #\N)
;;                              (:char v-char #\V)
;;                              (:char a-char #\A)
;;                              (:char u-char #\U))
;;           (with-work-areas ((m% :integer 1)
;;                             (isuppz% :integer (* 2 n)))
;;             (if vectors-p
;;                 (with-nv-output (z (* n n) z% real-lla-type)
;;                   (with-work-queries ((lwork work real-lla-type)
;;                                       (liwork iwork :integer))
;;                     (call-with-info-check procedure v-char a-char u-char n% a% n%
;;                                           null% null% null% null% ; VL,VU,IL,IU: not referenced
;;                                           abstol% m% 
;;                                           w% z% n% isuppz%
;;                                           work lwork iwork liwork info%))
;;                   (values w
;;                           (make-instance 'dense-matrix :nrow n :ncol n
;;                                          :data z)))
;;                 (progn
;;                   (with-work-queries ((lwork work real-lla-type)
;;                                       (liwork iwork :integer))
;;                     (call-with-info-check procedure n-char a-char u-char n% a% n%
;;                                           null% null% null% null% ; VL,VU,IL,IU: not referenced
;;                                           abstol% m% 
;;                                           w% null% n% isuppz%
;;                                           work lwork iwork liwork info%))
;;                   w))))))))

(defmethod eigen ((a hermitian-matrix) &key vectors-p check-real-p)
  ;; Uses simple driver, where "simple" means "silly collection of special cases".
  (declare (ignore check-real-p))       ; eigenvalues are always real
  (bind (((:slots-read-only (n nrow) (n2 ncol) (a-data data)) a)
         (lla-type (lla-type a-data))
         ((:values procedure real-lla-type complex-p) (lb-procedure-name2 'syev 'heev lla-type)))
    (assert (= n n2))
    (with-nv-output (w n w% real-lla-type) ; eigenvalues
      (with-fortran-atoms ((:integer n% n)
                           (:char nv-char (if vectors-p #\V #\N))
                           (:char u-char #\U)) ; upper triangle
        (if complex-p
            ;; *heev require an extra workspace argument
            (with-work-area (rwork% real-lla-type (max 1 (- (* 3 n) 2)))
              (if vectors-p
                  (with-nv-input-output (a-data v-data a% lla-type) ; eigenvectors end up here
                    (with-work-query (lwork% work% real-lla-type)
                      (call-with-info-check procedure nv-char u-char n% a% n% w% work% lwork% rwork% info%)
                      (values w
                              (make-instance 'dense-matrix :nrow n :ncol n
                                             :data v-data))))
                  (with-nv-input-copied (a-data a% lla-type) ; A is overwritten
                    (with-work-query (lwork% work% real-lla-type)
                      (call-with-info-check procedure nv-char u-char n% a% n% w% work% lwork% rwork% info%)
                      w))))
            (if vectors-p
                (with-nv-input-output (a-data v-data a% lla-type) ; eigenvectors end up here
                  (with-work-query (lwork% work% real-lla-type)
                    (call-with-info-check procedure nv-char u-char n% a% n% w% work% lwork% info%)
                    (values w
                            (make-instance 'dense-matrix :nrow n :ncol n
                                           :data v-data))))
                (with-nv-input-copied (a-data a% lla-type) ; A is overwritten
                  (with-work-query (lwork% work% real-lla-type)
                    (call-with-info-check procedure nv-char u-char n% a% n% w% work% lwork% info%)
                    w))))))))


;;;;
;;;; least squares calculations
;;;;

(defgeneric least-squares (a b)
  (:documentation "Return (values x qr ss nu), where x = argmin_x
L2norm( b-Ax ), solving a least squares problem, qr is the QR
decomposition of A, and SS is the sum of squares for each column of B.
B can have multiple columns, in which case x will have the same number
of columns, each corresponding to a different column of b.  nu is the
degrees of freedom."))

(defmethod least-squares ((a dense-matrix-like) (b dense-matrix-like))
  (set-restricted a)
  (set-restricted b)
  (bind (((:slots-read-only (m nrow) (n ncol) (a-data data)) a)
         ((:slots-read-only (m2 nrow) (nrhs ncol) (b-data data)) b)
         (common-type (smallest-common-target-type
                       (mapcar #'lla-type (list a-data b-data))))
         (procedure (lb-procedure-name 'gels common-type)))
    (assert (= m m2))
    (unless (<= n m)
      (error "A doesn't have enough columns for least squares"))
    (with-nv-input-output (a-data qr-data a% common-type) ; output: QR decomposition
      (with-nv-input-output (b-data x-data b% common-type) ; output: 
        (with-fortran-atoms ((:integer n% n)
                             (:integer m% m)
                             (:integer nrhs% nrhs)
                             (:char n-char #\N))
            (with-work-query (lwork% work% :double)
              (call-with-info-check procedure n-char m% n% nrhs% a% m% b% m%
                                    work% lwork% info%)))
        (values 
          (matrix-from-first-rows x-data m nrhs n)
          (make-instance 'qr :nrow m :ncol n
                         :data qr-data)
          (sum-last-rows x-data m nrhs n)
          (- m n))))))

;;; univariate versions of least squares: vector ~ vector, vector ~ matrix

(defmethod least-squares ((a dense-matrix-like) (b numeric-vector))
  (bind (((:values x qr ss nu) (least-squares a (vector->matrix-col b))))
    (values (matrix->vector x) qr (xref ss 0) nu)))

(defmethod least-squares ((a numeric-vector) (b numeric-vector))
  (bind (((:values x qr ss nu) (least-squares (vector->matrix-col a)
                                           (vector->matrix-col b))))
    (values (matrix->vector x) qr (xref ss 0) nu)))


;;;;
;;;; inverting matrices
;;;;

(defgeneric invert (a)
  (:documentation "Invert A."))

(defmethod invert ((a dense-matrix-like))
  (set-restricted a)
  (invert (lu a)))

(defmethod invert ((lu lu))
  (bind (((:slots-read-only (n nrow) (n2 ncol) ipiv (lu-data data)) lu)
	 (common-type (lla-type lu-data))
	 (procedure (lb-procedure-name 'getri common-type)))
    (assert (= n n2))
    (with-nv-input-output (lu-data inverse-data lu% common-type)
      (with-nv-input (ipiv ipiv% :integer)
        (with-fortran-atom (:integer n% n)
          (with-work-query (lwork% work% common-type)
            (call-with-info-check procedure n% lu% n% ipiv% work% lwork% info%))
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
	 (procedure (lb-procedure-name 'trtri common-type)))
    (assert (= n n2))
    (with-nv-input-output (a-data inv-data a% common-type)
      (with-fortran-atoms ((:char u-char (if upper-p #\U #\L))
                           (:char d-char (if unit-diag-p #\U #\N))
                           (:integer n% n))
        (call-with-info-check procedure u-char d-char n% a% n% info%))
      (make-instance result-class :nrow n :ncol n :data inv-data))))
  
(defmethod invert ((a upper-triangular-matrix))
  (invert-triangular% a t nil 'upper-triangular-matrix))

(defmethod invert ((a lower-triangular-matrix))
  (invert-triangular% a nil nil 'lower-triangular-matrix))

(defmethod invert ((a cholesky))
  (bind (((:slots-read-only (n nrow) (n2 ncol) (a-data data)) a)
	 (common-type (lla-type a-data))
	 (procedure (lb-procedure-name 'potri common-type)))
    (assert (= n n2))
    (with-nv-input-output (a-data inv-data a% common-type) ; output: upper triangular
      (with-fortran-atoms ((:integer n% n)
                           (:char u-char #\U))
        (call-with-info-check procedure u-char n% a% n% info%))
      (make-instance 'hermitian-matrix :nrow n :ncol n :data inv-data))))

;;;;
;;;;  utility functions for least squares
;;;;

(defun least-squares-raw-variance (qr)
  "Calculate the residual variance (basically, (X^T X)-1 ) from the qr
decomposition of X.  Return a CHOLESKY decomposition, which can be
used as a LOWER-TRIANGULAR-MATRIX for generating random draws. "
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
  (set-restricted a)
  (set-restricted b)
  (bind (((:slots-read-only (m nrow) (k ncol) (a-data data)) a)
         ((:slots-read-only (k2 nrow) (n ncol) (b-data data)) b)
         (common-type (smallest-common-target-type
                       (mapcar #'lla-type (list a-data b-data))))
         (lisp-type (lla-type->lisp-type common-type))
         (procedure (lb-procedure-name 'gemm common-type)))
    (assert (= k k2))
    (with-nv-input (a-data a% common-type)
      (with-nv-input (b-data b% common-type)
        (with-nv-output (c-data (* m n) c% common-type)
          (with-fortran-atoms ((:integer n% n)
                               (:integer k% k)
                               (:integer m% m)
                               (common-type alpha% (coerce alpha lisp-type))
                               (common-type z% (coerce 0 lisp-type))
                               (:char trans% #\N))
              (funcall procedure trans% trans% m% n% k% alpha% a% m% b% k%
                       z% c% m%))
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


(defmethod mm ((a numeric-vector) (b dense-matrix-like) &key (alpha 1))
  (matrix->vector (mm (vector->matrix-row a) b :alpha alpha)))

(defmethod mm ((a dense-matrix-like) (b numeric-vector) &key (alpha 1))
  (matrix->vector (mm a (vector->matrix-col b) :alpha alpha)))

(defmethod mm ((a numeric-vector) (b numeric-vector) &key (alpha 1))
  ;; a dot product, basically :-)
  (xref (mm (vector->matrix-row a) (vector->matrix-col b) :alpha alpha)
        0 0))
  

;;;;
;;;;  Updates for hermitian matrices
;;;;

;; (defun ntc-operation-char% (operation)
;;   "Return corresponding character for LAPACK operations."
;;   (ecase operation
;;     (:none #\N)
;;     (:transpose #\T)
;;     (:conjugate #\C)))

(defun mmx (a op-left-p &key (alpha 1) (beta 0) (c (lla-type a)))
  "Calculate alpha*A*op(A) + beta*C if OP-LEFT-P is NIL, and
  alpha*op(A)*A + beta*C otherwise.  ALPHA and BETA default to 1 and
  0, respectively.  op() is always conjugate transpose, but may be
  implemented as a transpose if A is real, in which case the two are
  equivalent.

  C can also be a symbol, denoting an LLA type, in which case a
  conformable matrix of that type will be created.  The default is the
  LLA-TYPE of A."
  (set-restricted a)
  (bind (((:slots-read-only (nrow-a nrow) (ncol-a ncol) (a-data data)) a)
         ((:values dim-c other-dim-a op-char)
          (if op-left-p
              (values ncol-a nrow-a #\C) ; C is good for real, too
              (values nrow-a ncol-a #\N))))
    ;; if C is a symbol, create the requested matrix
    (if (symbolp c)
        (setf c (make-matrix 'hermitian-matrix dim-c dim-c :lla-type c))
        (progn
          (check-type c hermitian-matrix)
          (assert (= dim-c (nrow c) (ncol c)))))
    (bind (((:slots-read-only (c-data data)) c)
           (common-type (smallest-common-target-type
                         (mapcar #'lla-type (list c-data c-data))))
           ((:values procedure real-lla-type) (lb-procedure-name2 'syrk 'herk common-type))
           (real-lisp-type (lla-type->lisp-type real-lla-type)))
      (with-nv-input (a-data a% common-type)
        (with-nv-input-output (c-data result-data c% common-type)
          (with-fortran-atoms ((:integer nrow-a)
                               (:integer dim-c)
                               (:integer other-dim-a)
                               (real-lla-type alpha (coerce alpha real-lisp-type))
                               (real-lla-type beta (coerce beta real-lisp-type))
                               (:char u-char #\U)
                               (:char op-char))
            (funcall procedure u-char op-char dim-c other-dim-a
                     alpha a% nrow-a beta c% dim-c))
          (make-instance 'hermitian-matrix :data result-data :nrow dim-c :ncol dim-c
                         :restricted-set-p nil))))))

;;;;
;;;;  Cholesky factorization
;;;;

(defgeneric cholesky (a)
  (:documentation "Cholesky factorization.  Only uses the lower
  triangle of a dense-matrix, and needs a PSD matrix."))

(defmethod cholesky ((a dense-matrix-like))
  (set-restricted a)
  (bind (((:slots-read-only (n nrow) (n2 ncol) data) a)
         (common-type (lla-type data))
         (procedure (lb-procedure-name 'potrf common-type)))
    (assert (= n n2))
    (with-nv-input-output (data cholesky-data data% common-type)
      (with-fortran-atoms ((:integer n% n)
                           (:char u-char #\L))
        (call-with-info-check procedure u-char n% data% n% info%))
      (make-instance 'cholesky :nrow n :ncol n
                     :data cholesky-data))))

(defmethod reconstruct ((mf cholesky))
  (mmx (take 'lower-triangular-matrix mf) nil))
