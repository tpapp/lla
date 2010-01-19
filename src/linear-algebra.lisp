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

                                                    
;;;;  matrix multiplication
;;; 
;;; The general matrix multiplication function is MMM.  It can be
;;; called with more than two arguments: it multiplies matrices from
;;; the left, using MM.  MM's first argument can also be a scalar.
;;;
;;; MM takes two matrices and an optional scalar.  Specialized methods
;;; may allow matrices to be replaced by symbols etc, in which case
;;; they have a special interpretation (eg multiply matrix by
;;; transpose).
;;;
;;; !! Various optimizations, deferred to the future:
;;; - compiler macros could notice (transpose matrix)
;;; - mmm could detect and collect scalars
;;; - methods for triangular matrices

(defun mmm (&rest matrices)
  (reduce #'mm matrices))

(defgeneric mm (a b &optional alpha)
  (:documentation "multiply A and B, also by the scalar
  alpha (defaults to 1)."))

(defmethod mm ((a numeric-vector) (b dense-matrix-like) &optional (alpha 1))
  (copy-nv (mm (vector->row a) b alpha)))

(defmethod mm ((a dense-matrix-like) (b numeric-vector) &optional (alpha 1))
  (copy-nv (mm a (vector->column b) alpha)))

(defmethod mm ((a numeric-vector) (b numeric-vector) &optional (alpha 1))
  ;; a dot product, basically
  (aref (elements (mm (vector->row a) (vector->column b) alpha)) 0))

(defmethod mm ((a dense-matrix-like) (b dense-matrix-like) &optional (alpha 1))
  (bind ((common-type (lb-target-type a b))
         (lisp-type (lla-type->lisp-type common-type))
         (procedure (lb-procedure-name 'gemm common-type)))
    (with-matrix-inputs (((a (m m%) (k k%)) a% common-type)
                         ((b k2 (n n%)) b% common-type))
      (assert (= k k2))
      (with-vector-output (c (* m n) c% common-type)
        (with-fortran-atoms ((common-type alpha% (coerce alpha lisp-type))
                             (common-type z% (coerce 0 lisp-type))
                             (:char trans% #\N))
          (funcall procedure trans% trans% m% n% k% alpha% a% m% b% k%
                   z% c% m%))
          (make-matrix* common-type m n c)))))

(defun mm-hermitian% (a op-left-p &optional (alpha 1))
  "Calculate alpha*A*op(A) if OP-LEFT-P is NIL, and alpha*op(A)*A
  otherwise.  ALPHA defaults to 1.  op() is always conjugate
  transpose, but may be implemented as a transpose if A is real, in
  which case the two are equivalent.  This function is meant to be
  used internally, and is *NOT EXPORTED*."
  (bind ((common-type (lla-type a))
         ((:values procedure real-type)
          (lb-procedure-name2 'syrk 'herk common-type)))
    (with-matrix-input ((a (nrow nrow%) ncol) a% common-type)
      (bind (((:values dim-c other-dim-a op-char)
              (if op-left-p
                  (values ncol nrow #\C) ; C is good for real, too
                  (values nrow ncol #\N)))
             (real-lisp-type (lla-type->lisp-type real-type)))
        (with-vector-output (c (expt dim-c 2) c% common-type)
          (with-fortran-atoms ((:integer dim-c% dim-c)
                               (:integer other-dim-a% other-dim-a)
                               (real-type alpha% (coerce alpha 
                                                         real-lisp-type))
                               (real-type beta% (coerce 0 real-lisp-type))
                               (:char u-char% #\U)
                               (:char op-char% op-char))
            (funcall procedure u-char% op-char% dim-c% other-dim-a%
                     alpha% a% nrow% beta% c% dim-c%))
          (make-matrix* common-type dim-c dim-c c :kind :hermitian))))))

(defmethod mm ((a dense-matrix-like) (b (eql t)) &optional (alpha 1))
  ;; A A^T
  (mm-hermitian% a nil alpha))

(defmethod mm ((a (eql t)) (b dense-matrix-like) &optional (alpha 1))
  ;; A^T A
  (mm-hermitian% b t alpha))

(defmethod mm ((a diagonal) (b dense-matrix-like) &optional (alpha 1))
  (bind ((common-type (common-target-type a b alpha))
         ((:slots-read-only (diagonal-elements elements)) a)
         ((:slots-read-only nrow ncol (matrix-elements elements)) b)
         (result (make-matrix nrow ncol common-type :kind (matrix-kind b)))
         ((:slots-read-only (result-elements elements)) result)
         (i 0))
    (assert (= (length diagonal-elements) nrow) () "Dimension mismatch.")
    (dotimes (col ncol)
      (dotimes (row nrow)
        (setf (aref result-elements i)
              (* (aref matrix-elements i) (aref diagonal-elements row) alpha))
        (incf i)))
    result))

(defmethod mm ((a dense-matrix-like) (b diagonal) &optional (alpha 1))
  (bind ((common-type (common-target-type a b alpha))
         ((:slots-read-only nrow ncol (matrix-elements elements)) a)
         ((:slots-read-only (diagonal-elements elements)) b)
         (result (make-matrix nrow ncol common-type :kind (matrix-kind a)))
         ((:slots-read-only (result-elements elements)) result)
         (i 0))
    (assert (= (length diagonal-elements) ncol) () "Dimension mismatch.")
    (dotimes (col ncol)
      (let ((d*alpha (* (aref diagonal-elements col) alpha)))
        (dotimes (row nrow)
          (setf (aref result-elements i)
                (* (aref matrix-elements i) d*alpha))
          (incf i))))
    result))


;;;; LU factorization

(defgeneric lu (a)
  (:documentation "LU decomposition of A"))

(defmethod lu ((a dense-matrix-like))
  (bind ((type (lla-type a))
         (procedure (lb-procedure-name 'getrf type)))
    (with-matrix-input ((a (m m%) (n n%) :output-to lu) a% type)
      (with-vector-output (ipiv (min m n) ipiv% :integer)
        (call-with-info-check procedure m% n% a% m% ipiv% info%)
        (make-instance 'lu 
                       :lu-matrix (make-matrix* type m n lu)
                       :ipiv (make-nv* :integer ipiv))))))

;;;; solving linear equations

(defgeneric solve (a b)
  (:documentation "Return X that solves AX=B."))

(defmethod solve (a (b numeric-vector))
  ;; Simply convert and call the method for a matrix.
  ;;
  ;; ??? check if the overhead is expensive, my hunch is that it
  ;; should be trivial -- Tamas
  (copy-nv (solve a (vector->column b))))

(defmethod solve ((a dense-matrix-like) (b dense-matrix-like))
  (bind ((common-type (lb-target-type a b))
         (procedure (lb-procedure-name 'gesv common-type)))
    (with-matrix-inputs (((a (n n%) n2 :copied) a% common-type)
                         ((b n3 (nrhs nrhs%) :output-to x-elements) b% common-type))
      (assert (= n n2 n3))
      (with-work-area (ipiv% :integer n) ; not saving IPIV
        (call-with-info-check procedure n% nrhs% a% n% ipiv% b% n% info%))
      (make-matrix* common-type n nrhs x-elements))))

(defmethod solve ((lu lu) (b dense-matrix-like))
  (bind (((:slots-read-only lu-matrix ipiv) lu)
         (common-type (lb-target-type lu-matrix b))
         (procedure (lb-procedure-name 'getrs common-type)))
    (with-matrix-inputs (((lu-matrix (n n%) n2) lu% common-type)
                         ((b n3 (nrhs nrhs%) :output-to x) b% common-type))
      (assert (= n n2 n3))
      (with-nv-input ((ipiv) ipiv% :integer)
        (with-fortran-atom (:char trans% #\N)
          (call-with-info-check procedure trans% n% nrhs% lu% n% ipiv% b% n% info%)))
      (make-matrix* common-type n nrhs x))))

(defgeneric invert (a &key &allow-other-keys)
  (:documentation "Invert A.  Usage note: inverting matrices is
  unnecessary and unwise in most cases, because it is numerically
  unstable.  If you are solving many Ax=b equations with the same A,
  use a matrix factorization like LU.  Most INVERT methods use a
  matrix factorization anyway."))

(defmethod invert ((a dense-matrix-like) &key)
  (invert (lu a)))

(defmethod invert ((lu lu) &key)
  (bind (((:slots-read-only lu-matrix ipiv) lu)
         (common-type (lla-type lu-matrix))
	 (procedure (lb-procedure-name 'getri common-type)))
    (with-matrix-input ((lu-matrix (n n%) n2 :output-to inverse) lu% common-type)
      (assert (= n n2))
      (with-nv-input ((ipiv) ipiv% :integer)
        (with-work-query (lwork% work% common-type)
          (call-with-info-check procedure n% lu% n% ipiv% work% lwork%
                                info%)))
      (make-matrix* common-type n n inverse))))

(defun invert-triangular% (a upper-p unit-diag-p result-kind)
  "Invert a dense (triangular) matrix using the LAPACK routine *TRTRI.
UPPER-P indicates if the matrix is in the upper or the lower triangle
of a (which needs to be a subtype of dense-matrix-like, but the type
information is not used), UNIT-DIAG-P indicates whether the diagonal
is supposed to consist of 1s.  *For internal use, NOT EXPORTED*."
  (bind ((common-type (lla-type a))
	 (procedure (lb-procedure-name 'trtri common-type)))
    (with-matrix-input ((a (n n%) n2 :output-to inverse) a% common-type)
      (assert (= n n2))
      (with-fortran-atoms ((:char u-char% (if upper-p #\U #\L))
                           (:char d-char% (if unit-diag-p #\U #\N)))
        (call-with-info-check procedure u-char% d-char% n% a% n% info%))
      (make-matrix* common-type n n inverse :kind result-kind))))
  
(defmethod invert ((a upper-triangular-matrix) &key)
  (invert-triangular% a t nil :upper-triangular))

(defmethod invert ((a lower-triangular-matrix) &key)
  (invert-triangular% a nil nil :lower-triangular))

(defmethod invert ((cholesky cholesky) &key)
  (bind (((:slots-read-only r-matrix) cholesky)
         (common-type (lla-type r-matrix))
	 (procedure (lb-procedure-name 'potri common-type)))
    (with-matrix-input ((r-matrix (n n%) n2 :output-to inverse) r% common-type)
      (assert (= n n2))
      (with-fortran-atom (:char u-char% #\U)
        (call-with-info-check procedure u-char% n% r% n% info%))
      (make-matrix* common-type n n inverse :kind :hermitian))))

(defmethod invert ((d diagonal) &key (tolerance 0))
  "For pseudoinverse, suppressing diagonal elements below TOLERANCE
\(if given, otherwise / is just used without any checking."
  (bind (((:values elements type) (float-elements% d))
         (zero (coerce* 0 type)))
    (make-diagonal* type (map (nv-array-type type)
                              (cond
                                ((null tolerance) #'/)
                                ((and (numberp tolerance) (<= 0 tolerance))
                                 (let ((tolerance (coerce* tolerance
                                                           type)))
                                   (lambda (x)
                                     (if (<= (abs x) tolerance) zero (/ x)))))
                                ((and (numberp tolerance) (zerop tolerance))
                                 (lambda (x)
                                   (if (zerop x) zero (/ x))))
                                (t (error "Invalid tolerance argument.")))
                              elements))))

(defgeneric eigen (a &key vectors-p check-real-p &allow-other-keys)
  (:documentation "Calculate the eigenvalues and optionally the right
eigenvectors of a matrix (as columns).  Return (values eigenvalues
eigenvectors).  If check-real-p, eigenvalues of real matrices are
checked for an imaginary part and returned with the appropriate
type (compex or not).  Complex conjugate pairs of eigenvalues appear
consecutively with the eigenvalue having the positive imaginary part
first."))

(defun eigen-dense-real% (a vectors-p check-real-p)
  "Eigenvalues and vectors for dense, real matrices."
  ;; This is probably the hairiest function in the whole library.
  ;; S/DEEV returns real and imaginary parts for eigenvectors and
  ;; eigenvalues separately, so they have to be "zipped" together,
  ;; with two utility functions.
  (bind (((:values real-type complex-type)
          (ecase (lla-type a)
            ((:single :integer) (values :single :complex-single))
            (:double (values :double :complex-double))))
         (procedure (lb-procedure-name 'geev real-type)))
    (with-matrix-input ((a (n n%) n2 :copied) a% real-type) ; overwritten
      (assert (= n n2))
      (with-work-area (w% real-type (* 2 n)) ; eigenvalues, will be zipped
        (let (;; imaginary part
              (wi% (inc-pointer w% (* n (foreign-type-size real-type)))))
          (with-fortran-atoms ((:char n-char% #\N)
                               (:char v-char% #\V))
            (if vectors-p
                (with-vector-output (vr (expt n 2) vr% :double)
                  (with-work-query (lwork work :double)
                    (call-with-info-check procedure n-char% v-char% n% a% n% 
                                          w% wi% ; eigenvalues
                                          (null-pointer) n% vr% n% ; eigenvectors
                                          work lwork info%))
                  (bind (((:values eigenvalues complex-p)
                          (zip-eigenvalues w% n real-type complex-type check-real-p)))
                    (values eigenvalues
                            (if complex-p
                                (make-matrix* complex-type n n
                                              (zip-eigenvectors w% vr% n real-type
                                                                complex-type))
                                (make-matrix* real-type n n vr)))))
                (with-work-query (lwork work :double)
                  (call-with-info-check procedure n-char% n-char% n% a% n% 
                                        w% wi%   ; eigenvalues
                                        (null-pointer) n% (null-pointer) n% ; eigenvectors
                                        work lwork info%)
                  (zip-eigenvalues w% n real-type complex-type check-real-p)))))))))

(defun eigen-dense-complex% (a vectors-p)
  "Eigenvalues and vectors for dense, complex matrices."
  (bind ((complex-type (lla-type a))
         (procedure (lb-procedure-name 'geev complex-type)))
    (assert (lla-complex-p complex-type) () "This function only handles complex matrices.")
    (with-matrix-input ((a (n n%) n2 :copied) a% complex-type)
      (assert (= n n2))
      (with-vector-output (w n w% complex-type)
        (with-work-area (rwork complex-type n)
          (with-fortran-atoms ((:char n-char% #\N)
                               (:char v-char% #\V))
            (if vectors-p
                (with-vector-output (vr (expt n 2) vr% complex-type)
                  (with-work-query (lwork work complex-type)
                    (call-with-info-check procedure n-char% v-char% n% a% n% w%
                                          (null-pointer) n% vr% n% work lwork rwork info%)
                    (values (make-nv* complex-type w)
                            (make-matrix* complex-type n n vr))))
                (with-work-query (lwork work complex-type)
                  (call-with-info-check procedure n-char% n-char% n% a% n% w%
                                        (null-pointer) n% (null-pointer) n% 
                                        work lwork rwork info%)
                  (make-nv* complex-type w)))))))))

(defmethod eigen ((a dense-matrix-like) &key vectors-p check-real-p)
  ;; Unfortunately, the LAPACK interface for real and complex cases is
  ;; different.
  (case (lla-type a)
    ((:integer :single :double)
       (eigen-dense-real% a vectors-p check-real-p))
    ((:complex-single :complex-double)
       (eigen-dense-complex% a vectors-p))))

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
  ;; Uses simple driver, where "simple" means "silly collection of
  ;; special cases for all combinations".
  (declare (ignore check-real-p))       ; eigenvalues are always real
  (bind ((common-type (lla-type a))
         ((:values procedure real-type complex-p)
          (lb-procedure-name2 'syev 'heev common-type))
         (n (nrow a))
         ((:values w v)
          (with-fortran-atoms ((:char nv-char% (if vectors-p #\V #\N))
                               (:char u-char% #\U)) ; upper triangle
            (if complex-p
                 ;; *heev require an extra workspace argument
                (if vectors-p
                    ;; complex, with vectors
                    (with-matrix-input ((a (n n%) n2 :output-to v) a% common-type nil)
                      (assert (= n n2))
                      (with-vector-output (w n w% real-type) ; eigenvalues
                        (with-work-area (rwork% real-type (max 1 (- (* 3 n) 2)))
                          (with-work-query (lwork% work% real-type)
                            (call-with-info-check procedure nv-char% u-char% n% a% n% w%
                                                  work% lwork% rwork% info%))
                          (values w v))))
                    ;; complex, without vectors
                    (with-matrix-input ((a (n n%) n2 :copied) a% common-type nil)
                      (assert (= n n2))
                      (with-vector-output (w n w% real-type) ; eigenvalues
                        (with-work-area (rwork% real-type (max 1 (- (* 3 n) 2)))
                          (with-work-query (lwork% work% real-type)
                            (call-with-info-check procedure nv-char% u-char% n% a% n% w%
                                                  work% lwork% rwork% info%))
                          (values w nil)))))
                (if vectors-p
                    ;; real, with vectors
                    (with-matrix-input ((a (n n%) n2 :output-to v) a% common-type nil)
                      (assert (= n n2))
                      (with-vector-output (w n w% real-type) ; eigenvalues
                        (with-work-query (lwork% work% real-type)
                          (call-with-info-check procedure nv-char% u-char% n% a% n% w%
                                                work% lwork% info%))
                        (values w v)))
                    ;; real, without vectors
                    (with-matrix-input ((a (n n%) n2 :output-to v) a% common-type nil)
                      (assert (= n n2))
                      (with-vector-output (w n w% real-type) ; eigenvalues
                        (with-work-query (lwork% work% real-type)
                          (call-with-info-check procedure nv-char% u-char% n% a% n% w%
                                                work% lwork% info%))
                        (values w nil))))))))
    (values (make-nv* real-type w)
            (if v
                (make-matrix* real-type n n v)
                nil))))

;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; moving barrier -- things below are to be rewritten
;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;
;;;; least squares calculations
;;;;

(defgeneric least-squares (a b)
  (:documentation "Return (values x qr ss nu), where x = argmin_x
L2norm( b-Ax ), solving a least squares problem, QR is the QR
decomposition of A, and SS is the sum of squares for each column of B.
B can have multiple columns, in which case X will have the same number
of columns, each corresponding to a different column of B.  nu is the
degrees of freedom."))

(defmethod least-squares ((a dense-matrix-like) (b dense-matrix-like))
  (bind ((common-type (lb-target-type a b))
         (procedure (lb-procedure-name 'gels common-type)))
    (with-matrix-inputs (((a (m m%) (n n%) :output-to qr) a% common-type)
                         ((b m2 (nrhs nrhs%) :output-to x) b% common-type))
      (assert (= m m2))
      (unless (<= n m)
        (error "A doesn't have enough columns for least squares"))
      (with-fortran-atoms ((:char n-char% #\N))
        (with-work-query (lwork% work% :double)
          (call-with-info-check procedure n-char% m% n% nrhs% a% m% b% m%
                                work% lwork% info%)))
      (values 
        (matrix-from-first-rows x common-type m nrhs n)
        (make-instance 'qr :qr-matrix (make-matrix* common-type m n qr))
        (sum-last-rows x common-type m nrhs n)
        (- m n)))))


;;; univariate versions of least squares: vector ~ vector, vector ~ matrix

(defmethod least-squares ((a dense-matrix-like) (b numeric-vector))
  (bind (((:values x qr ss nu) (least-squares a (vector->column b))))
    (values (copy-nv x) qr (aref (elements ss) 0) nu)))

(defmethod least-squares ((a numeric-vector) (b numeric-vector))
  (bind (((:values x qr ss nu) (least-squares (vector->column a)
                                              (vector->column b))))
    (values (aref (elements x) 0) qr (aref (elements ss) 0) nu)))

(defun least-squares-xxinverse (qr)
  "Calculate (X^T X)-1 (which is used for calculating the variance of
estimates) from the qr decomposition of X.  Return a CHOLESKY
decomposition.  Note: the R-MATRIX of the cholesky decomposition can
be used to generate random draws, etc."
  ;; Notes: X = QR, thus X^T X = R^T Q^T Q R = R^T R because Q is
  ;; orthogonal.  Then we do as if calculating the inverse of a matrix
  ;; using its Cholesky factorization.
  (with-slots (nrow ncol) (qr-matrix qr)
    (assert (<= ncol nrow))
    (invert (make-instance 'cholesky :r-matrix (factorization-component qr :R)))))

;;;; Cholesky factorization

(defgeneric cholesky (a)
  (:documentation "Cholesky factorization.  Only uses the lower
  triangle of a dense-matrix, and needs a PSD hermitian matrix."))

(defmethod cholesky ((a hermitian-matrix))
  (bind ((common-type (lla-type a))
         (procedure (lb-procedure-name 'potrf common-type)))
    (with-matrix-input ((a (n n%) n2 :output-to cholesky) a% common-type)
      (assert (= n n2))
      (with-fortran-atoms ((:char u-char% #\L))
        (call-with-info-check procedure u-char% n% a% n% info%))
      (make-instance 'cholesky :r-matrix (make-matrix* common-type n n cholesky 
                                                       :kind :upper-triangular)))))

(defmethod reconstruct ((mf cholesky))
  (mm t (r-matrix mf)))



;;; SVD

(defun svd (a left-vector-spec right-vector-spec)
  "Singular value decomposition.  Valid vector specifications
are :NONE, :SINGULAR (singular vectors only) and :ALL.  Return values
S (singular values, descending order, as a DIAGONAL), U (left singular
vectors, DENSE-MATRIX), VT ([conjugate] transpose of right singular
vectors, DENSE-MATRIX)."
  (bind ((type (lla-type a))
         ((:values procedure real-type complex-p) (lb-procedure-name2 'gesvd 'gesvd type)))
    (with-matrix-input ((a (m m%) (n n%) :copied) a% type)
      (bind ((min-mn (min m n))
             ((:values u-ncol jobu) (ecase left-vector-spec
                                      (:none (values 0 #\N))
                                      (:singular (values min-mn #\S))
                                      (:all (values m #\A))))
             ((:values vt-nrow jobvt) (ecase right-vector-spec
                                        (:none (values 0 #\N))
                                        (:singular (values min-mn #\S))
                                        (:all (values n #\A)))))
        (with-fortran-atoms ((:char jobu% jobu)
                             (:char jobvt% jobvt)
                             (:integer vt-nrow% vt-nrow))
          (with-vector-output (s min-mn s% real-type)
            (bind (((:values u vt)
                    (with-vector-outputs ((u (* m u-ncol ) u% type)
                                          (vt (* vt-nrow n) vt% type))
                      (with-work-query (lwork% work% type)
                        (if complex-p
                            (with-work-area (rwork% real-type (* 5 min-mn))
                              (call-with-info-check procedure jobu% jobvt% m% n% a% m% s%
                                                    u% m% vt% vt-nrow% work% lwork% rwork% info%))
                            (call-with-info-check procedure jobu% jobvt% m% n% a% m% s%
                                                  u% m% vt% vt-nrow% work% lwork% info%)))
                      (values u vt))))
              (values (make-diagonal* type s)
                      (if (zerop u-ncol)
                          nil
                          (make-matrix* type m u-ncol u))
                      (if (zerop vt-nrow)
                          nil
                          (make-matrix* type vt-nrow n vt))))))))))
