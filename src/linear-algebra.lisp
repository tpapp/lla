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

(defun mm-nv-row% (a b alpha)
  "Matrix multiplication, with A converted to a row matrix, and the
product converted back to a numeric-vector."
  (copy-nv (mm (vector->row a) b alpha)))

(defun mm-nv-column% (a b alpha)
  "Matrix multiplication, with B converted to a column matrix, and the
product converted back to a numeric-vector."
  (copy-nv (mm a (vector->column b) alpha)))

(defmethod mm ((a numeric-vector) (b dense-matrix-like) &optional (alpha 1))
  (mm-nv-row% a b alpha))

(defmethod mm ((a dense-matrix-like) (b numeric-vector) &optional (alpha 1))
  (mm-nv-column% a b alpha))

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
         (result (make-matrix common-type nrow ncol :kind (matrix-kind b)))
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
         (result (make-matrix common-type nrow ncol :kind (matrix-kind a)))
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

(defmethod mm ((a numeric-vector) (b diagonal) &optional (alpha 1))
  (mm-nv-row% a b alpha))

(defmethod mm ((a diagonal) (b numeric-vector) &optional (alpha 1))
  (mm-nv-column% a b alpha))


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

;;;; Hermitian factorization

(defun hermitian-factorization (a &key (component :U))
  (check-type a hermitian-matrix)
  (bind ((type (lla-type a))
         (procedure (lb-procedure-name2 'sytrf 'hetrf type))
         ((:values u-char kind set-restricted-p) (ecase component
                                                   (:U (values #\U :upper-triangular nil))
                                                   (:D (values #\L :lower-triangular t)))))
    (with-matrix-input ((a (n n%) n2 :output-to factor) a% type set-restricted-p)
      (assert (= n n2) () "Hermitian matrix is not square.")
      (with-vector-output (ipiv n ipiv% :integer)
        (with-fortran-atom (:char u-char% u-char)
          (with-work-query (lwork% work% type)
            (call-with-info-check procedure u-char% n% a% n% ipiv% work% lwork% info%)))
        (make-instance 'hermitian-factorization 
                       :factor (make-matrix* type n n factor :kind kind)
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

(defmethod solve ((cholesky cholesky) (b dense-matrix-like))
  (bind (((:slots-read-only factor) cholesky)
         (common-type (lla-type factor))
	 (procedure (lb-procedure-name 'potrs common-type)))
    (with-matrix-inputs (((factor (n n%) n2) factor% common-type)
                         ((b n3 (nrhs nrhs%) :output-to x) b% common-type))
      (assert (= n n2 n3))
      (with-fortran-atom (:char u-char% (etypecase factor
                                          (upper-triangular-matrix #\U)
                                          (lower-triangular-matrix #\L)))
        (call-with-info-check procedure u-char% n% nrhs% factor% n% b% n% info%))
      (make-matrix* common-type n nrhs x))))

(defmethod solve ((hermitian-factorization hermitian-factorization) (b dense-matrix-like))
  (bind (((:slots-read-only factor ipiv) hermitian-factorization)
         (common-type (lla-type factor))
	 (procedure (lb-procedure-name2 'sytrs 'hetrs common-type)))
    (with-matrix-inputs (((factor (n n%) n2) factor% common-type)
                         ((b n3 (nrhs nrhs%) :output-to x) b% common-type))
      (assert (= n n2 n3))
      (with-nv-input ((ipiv) ipiv% :integer)
        (with-fortran-atom (:char u-char% (etypecase factor
                                            (upper-triangular-matrix #\U)
                                            (lower-triangular-matrix #\L)))
          (call-with-info-check procedure u-char% n% nrhs% factor% n% ipiv% b% n% info%))
        (make-matrix* common-type n n x)))))

(defun trsm% (a b side transpose-a? &optional (alpha 1))
  "Wrapper for BLAS routine xTRSM.  Calculates op(A^-1) B (if SIDE
is :LEFT) or B op(A^-1) (if SIDE is :RIGHT).  A has to be a triangular
matrix.  transpose-a? determines whether op(A) is A^T or A.  The
result is multiplied by ALPHA."
  (bind ((common-type (lb-target-type a b))
         (procedure (lb-procedure-name 'trsm common-type)))
    (with-matrix-inputs (((a (a-n lda%) a-m) a% common-type)
                         ((b (m m%) (n n%) :output-to result) b% common-type))
      (bind ((side (ecase side
                     (:left
                        (assert (= a-m m))
                        #\L)
                     (:right
                        (assert (= a-n n))
                        #\R))))
        (with-fortran-atoms ((:char side% side)
                             (:char uplo% (ecase (matrix-kind a)
                                            (:lower-triangular #\L)
                                            (:upper-triangular #\U)))
                             (:char transa% (if transpose-a?
                                                #\C
                                                #\N))
                             (:char diag% #\N)
                             (common-type alpha% (coerce* alpha common-type)))
          (funcall procedure side% uplo% transa% diag%
                   m% n% alpha% a% lda% b% m%)
          (make-matrix* common-type m n result))))))

(defmethod solve ((a lower-triangular-matrix) b)
  (trsm% a b :left nil))

(defmethod solve ((a upper-triangular-matrix) b)
  (trsm% a b :left nil))

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

(defmethod invert ((a hermitian-matrix) &key)
   (invert (hermitian-factorization a)))

(defmethod invert ((hf hermitian-factorization) &key)
  ;; If the FACTOR is lower triangular, we need to transpose it, as
  ;; hermitian matrices always store the upper triangle.
  (bind (((:slots-read-only factor ipiv) hf)
         (factor (aetypecase factor
                   (lower-triangular-matrix (transpose it))
                   (upper-triangular-matrix it)))
         (common-type (lla-type factor))
         (procedure (lb-procedure-name2 'sytri 'hetri common-type)))
    (with-matrix-input ((factor (n n%) n2 :output-to inverse) factor% common-type)
      (assert (= n n2))
      (with-work-area (work% common-type n)
        (with-nv-input ((ipiv) ipiv% :integer)
          (with-fortran-atom (:char u-char% #\U)
            (call-with-info-check procedure u-char% n% factor% n% ipiv% work% info%))))
      (make-matrix* common-type n n inverse :kind :hermitian))))

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
  ;; If the FACTOR of CHOLESKY is lower triangular, we need to
  ;; transpose it, as hermitian matrices always store the upper
  ;; triangle.
  (bind ((factor (aetypecase (factor cholesky)
                   (lower-triangular-matrix (transpose it))
                   (upper-triangular-matrix it)))
         (common-type (lla-type factor))
         (procedure (lb-procedure-name 'potri common-type)))
    (with-matrix-input ((factor (n n%) n2 :output-to inverse) factor% common-type)
      (assert (= n n2))
      (with-fortran-atom (:char u-char% #\U)
        (call-with-info-check procedure u-char% n% factor% n% info%))
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
  (declare (optimize debug))
  (bind (((:values real-type complex-type)
          (ecase (lla-type a)
            ((:single :integer) (values :single :complex-single))
            (:double (values :double :complex-double))))
         (procedure (lb-procedure-name 'geev real-type)))
    (with-matrix-input ((a (n n%) n2 :copied) a% real-type) ; overwritten
      (assert (= n n2))
      (with-work-area (w% real-type (* 2 n)) ; eigenvalues, will be zipped
        (let (;; imaginary part
              (wi% (inc-pointer w% (* n (foreign-size* real-type)))))
          (with-fortran-atoms ((:char n-char% #\N)
                               (:char v-char% #\V))
            (if vectors-p
                (with-vector-output (vr (expt n 2) vr% real-type)
                  (with-work-query (lwork work real-type)
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
                (with-work-query (lwork% work% real-type)
                  (call-with-info-check procedure n-char% n-char% n% a% n% 
                                        w% wi%   ; eigenvalues
                                        (null-pointer) n% (null-pointer) n% ; eigenvectors
                                        work% lwork% info%)
                  (zip-eigenvalues w% n real-type complex-type check-real-p)))))))))

(defun eigen-dense-complex% (a vectors-p)
  "Eigenvalues and vectors for dense, complex matrices."
  (bind ((complex-type (lla-type a))
         (procedure (lb-procedure-name 'geev complex-type)))
    (assert (lla-complex-p complex-type) () "This function only handles complex matrices.")
    (with-matrix-input ((a (n n%) n2 :copied) a% complex-type)
      (assert (= n n2))
      (with-vector-output (w n w% complex-type)
        (with-work-area (rwork% complex-type n)
          (with-fortran-atoms ((:char n-char% #\N)
                               (:char v-char% #\V))
            (if vectors-p
                (with-vector-output (vr (expt n 2) vr% complex-type)
                  (with-work-query (lwork% work% complex-type)
                    (call-with-info-check procedure n-char% v-char% n% a% n% w%
                                          (null-pointer) n% vr% n% work% lwork% rwork% info%))
                  (values (make-nv* complex-type w)
                          (make-matrix* complex-type n n vr)))
                (with-work-query (lwork% work% complex-type)
                  (call-with-info-check procedure n-char% n-char% n% a% n% w%
                                        (null-pointer) n% (null-pointer) n% 
                                        work% lwork% rwork% info%)
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
;;;  !!!! interface has changed quite a bit, clean up code below if
;;;  !!!! ever used
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

;;;;
;;;; least squares calculations
;;;;

(defgeneric least-squares (y x)
  (:documentation "Return (values b qr ss nu), where beta = argmin_b
L2norm( y-Xb ), solving a least squares problem, QR is the QR
decomposition of A, and SS is the sum of squares for each column of Y.
Y can have multiple columns, in which case X will have the same number
of columns, each corresponding to a different column of Y.  NU is the
degrees of freedom."))

(defmethod least-squares ((y dense-matrix-like) (x dense-matrix-like))
  ;; Note: the naming convention (y,X,b) is different from LAPACK's
  ;; (b,A,x).  Sorry if this creates confusion, I decided to follow
  ;; standard statistical notation.
  (bind ((common-type (lb-target-type x y))
         (procedure (lb-procedure-name 'gels common-type)))
    (with-matrix-inputs (((x (m m%) (n n%) :output-to qr) x% common-type)
                         ((y m2 (nrhs nrhs%) :output-to b) y% common-type))
      (assert (= m m2))
      (unless (<= n m)
        (error "A doesn't have enough columns for least squares"))
      (with-fortran-atoms ((:char n-char% #\N))
        (with-work-query (lwork% work% :double)
          (call-with-info-check procedure n-char% m% n% nrhs% x% m% y% m%
                                work% lwork% info%)))
      (values 
        (matrix-from-first-rows common-type b n nrhs m)
        (make-instance 'qr :qr-matrix (make-matrix* common-type m n qr))
        (sum-last-rows common-type b m nrhs n)
        (- m n)))))


;;; univariate versions of least squares: vector ~ vector, vector ~ matrix

(defmethod least-squares ((y numeric-vector) (x dense-matrix-like))
  (bind (((:values b qr ss nu) (least-squares (vector->column y) x)))
    (values (copy-nv b) qr (aref (elements ss) 0) nu)))

(defmethod least-squares ((y numeric-vector) (x numeric-vector))
  (bind (((:values b qr ss nu) (least-squares (vector->column y)
                                              (vector->column x))))
    (values (aref (elements b) 0) qr (aref (elements ss) 0) nu)))

(defun least-squares-xx-inverse (qr)
  "Calculate (X^T X)-1 (which is used for calculating the variance of
estimates) from the qr decomposition of X.  Return a CHOLESKY
decomposition.  Note: the FACTOR of the cholesky decomposition can be
used to generate random draws, etc."
  ;; Notes: X = QR, thus X^T X = R^T Q^T Q R = R^T R because Q is
  ;; orthogonal.  Then we do as if calculating the inverse of a matrix
  ;; using its Cholesky factorization.
  (with-slots (nrow ncol) (qr-matrix qr)
    (assert (<= ncol nrow))
    (invert (make-instance 'cholesky :factor (component qr :R)))))


;;;; constrained-least-squares

(defgeneric constrained-least-squares (y x z w)
  (:documentation "Solve the (linearly) constrained least squares
  problem min_b |y-Xb|_2 subject to Zx=w.  Return b."))

(defmethod constrained-least-squares ((y numeric-vector) (x dense-matrix-like)
                                      (z dense-matrix-like) (w numeric-vector))
  "Solve the (linearly) constrained least squares problem min_b |y-Xb|_2
  subject to Zx=w."
  ;; Note: mapping between the function parameters/variables and
  ;; LAPACK counterparts is as follows: y->c, X->A, x->b, Z->B, w->d,
  (bind ((common-type (lb-target-type y x z w))
         (procedure (lb-procedure-name 'gglse common-type)))
    ;; !! after call, x and z contain decompositions, currently not collected
    (with-matrix-inputs (((x (m m%) (n n%) :copied) x% common-type)
                         ((z (p p%) n2 :copied) z% common-type))
      (assert (= n n2) () "Dimension mismatch between z and x")
      (assert (= m (xsize y)) () "Dimension mismatch between x and y")
      (assert (= p (xsize w)) () "Dimension mismatch between z and w")
      ;; !! after call, y contains sum of squares, currently not collected
      (with-nv-inputs (((y :copied) y% common-type)
                       ((w :copied) w% common-type))
        (with-vector-output (b n b% common-type)
          (with-work-query (lwork% work% common-type)
            (call-with-info-check procedure m% n% p% x% m% z% p% y% w% b% work% lwork% info%))
          (make-nv* common-type b))))))


;;;; Cholesky factorization

(defgeneric cholesky (a &optional component)
  (:documentation "Cholesky factorization.  Component is :L or :U (default)."))

(defmethod cholesky ((a hermitian-matrix) &optional (component :U))
  (bind ((common-type (lla-type a))
         (procedure (lb-procedure-name 'potrf common-type))
         ((:values u-char kind) (ecase component
                                  (:U (values #\U :upper-triangular))
                                  (:L (values #\L :lower-triangular)))))
    (with-matrix-input ((a (n n%) n2 :output-to cholesky) a% common-type)
      (assert (= n n2))
      (with-fortran-atoms ((:char u-char% u-char))
        (call-with-info-check procedure u-char% n% a% n% info%))
      (make-instance 'cholesky :factor (make-matrix* common-type n n cholesky 
                                                     :kind kind)))))

(defmethod reconstruct ((mf cholesky))
  (bind (((:slots-read-only factor) mf))
    (etypecase factor
      (upper-triangular-matrix (mm t factor))
      (lower-triangular-matrix (mm factor t)))))



;;; SVD

(defun svd (a &key (left :none) (right :none))
  "Singular value decomposition.  Valid vector specifications
are :NONE, :SINGULAR (singular vectors only) and :ALL.  Return values
S (singular values, descending order, as a DIAGONAL), U (left singular
vectors, DENSE-MATRIX), VT ([conjugate] transpose of right singular
vectors, DENSE-MATRIX)."
  (bind ((type (lla-type a))
         ((:values procedure real-type complex-p)
          (lb-procedure-name2 'gesvd 'gesvd type)))
    (with-matrix-input ((a (m m%) (n n%) :copied) a% type)
      (bind ((min-mn (min m n))
             ((:values u-ncol jobu) (ecase left
                                      (:none (values 1 #\N)) ; LAPACK needs >=1
                                      (:singular (values min-mn #\S))
                                      (:all (values m #\A))))
             ((:values vt-nrow jobvt) (ecase right
                                        (:none (values 1 #\N)) ; LAPACK needs >=1
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
                        (with-lapack-traps-masked 
                          (if complex-p
                              (with-work-area (rwork% real-type (* 5 min-mn))
                                (call-with-info-check 
                                 procedure jobu% jobvt% m% n% a% m% s%
                                 u% m% vt% vt-nrow% work% lwork% rwork% info%))
                              (call-with-info-check 
                               procedure jobu% jobvt% m% n% a% m% s%
                               u% m% vt% vt-nrow% work% lwork% info%))))
                      (values u vt))))
              (values (make-diagonal* type s)
                      (if (eq left :none)
                          nil
                          (make-matrix* type m u-ncol u))
                      (if (eq right :none)
                          nil
                          (make-matrix* type vt-nrow n vt))))))))))

;;; trace

(defgeneric tr (a)
  (:documentation "Trace of a square matrix.")
  (:method ((a dense-matrix-like))
    (check-type a square-matrix)
    (bind (((:lla-matrix a) a))
      (iter
        (for i :from 0 :below (nrow a))
        (summing (a (a-index i i))))))
  (:method ((a diagonal))
    (reduce #'+ (elements a))))

;; rank

(defun rank (matrix &key threshold (logrc-threshold 0.5))
  "Calculate the rank of a matrix using singular values - only those
above THRESHOLD in absolute value are counted.  If NIL, THRESHOLD is
determined automatically as (* (MAX NROW NCOL) FLOAT-EPSILON) times
the largest absolute singular value.  If (ABS (LOG (/ NROW NCOL))) is
above LOGRC-THRESHOLD, use (MM MATRIX T) or its transpose for faster
calculations; this can be disabled by setting LOG-THRESHOLD to NIL.
Returns the absolute values of the singular values as the second
value."
  (check-type logrc-threshold (or null (and number (satisfies plusp))))
  (bind (((:accessors-r/o nrow ncol) matrix)
         (ratio (log (/ nrow ncol)))
         ((:values matrix squared-p)
          (cond 
            (;; fewer rows than columns
             (aand logrc-threshold (< ratio (- it)))
             (values (mm matrix t) t))
            (;; more rows than columns
             (aand logrc-threshold (< it ratio))
             (values (mm t matrix) t))
            (;; no matrix multiplication
             t
             (values matrix nil))))
         (d (svd matrix))
         (d-elements (elements d)))
    (map-into d-elements (if squared-p
                             (lambda (x)
                               (sqrt (abs x)))
                             #'abs)
              d-elements)
    (let ((threshold (aif threshold
                          it
                          (* (max nrow ncol) (epsilon* (lla-type d))
                             (reduce #'max d-elements)))))
    (values 
      (count-if (lambda (x) (<= threshold (abs x)))
                (elements d))
      d))))
