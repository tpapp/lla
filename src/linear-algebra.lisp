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
;;; the left, using MM.
;;;
;;; MM takes two matrices and an optional scalar.  Specialized methods
;;; may allow matrices to be replaced by symbols etc, in which case
;;; they have a special interpretation (eg T stands for "multiply
;;; matrix by conjugate transpose").
;;;
;;; !! Various optimizations, deferred to the future:
;;; - compiler macros could notice (transpose matrix)
;;; - mmm could detect and collect scalars
;;; - methods for triangular matrices

(defgeneric mm (a b &optional alpha)
  (:documentation "multiply A and B, also by the scalar
  alpha (defaults to 1)."))

(defun mmm (&rest matrices)
  (reduce #'mm matrices))

(defun mm-nv-row% (a b alpha)
  "Matrix multiplication, with A converted to a row matrix, and the
product converted back to a numeric-vector."
  (elements (mm (as-row a) b alpha)))

(defun mm-nv-column% (a b alpha)
  "Matrix multiplication, with B converted to a column matrix, and the
product converted back to a numeric-vector."
  (elements (mm a (as-column b) alpha)))

(defmethod mm ((a vector) (b dense-matrix-like) &optional (alpha 1))
  (mm-nv-row% a b alpha))

(defmethod mm ((a dense-matrix-like) (b vector) &optional (alpha 1))
  (mm-nv-column% a b alpha))

;;; for numeric vectors, mm is the dot product; !!! need to do this
;;; decently one day using BLAS

(defmethod mm ((a vector) (b vector) &optional (alpha 1))
  (as-scalar% (elements (mm (as-row a) (as-column b) alpha))))

(defmethod mm ((a vector) (b (eql t)) &optional (alpha 1))
  (mm a a alpha))

(defmethod mm ((a (eql t)) (b vector) &optional (alpha 1))
  (mm b b alpha))

(defmethod mm ((a dense-matrix-like) (b dense-matrix-like)
               &optional (alpha 1))
  (lb-call ((common-type (lb-target-type a b alpha))
            (procedure (lb-procedure-name common-type gemm))
            ((:matrix a% common-type (m m%) (k k%)) a)
            ((:matrix b% common-type k2 (n n%)) b)
            ((:output c% common-type c) (* m n))
            ((:atom alpha% common-type) (coerce* alpha common-type))
            ((:atom z% common-type) (zero* common-type))
            ((:char trans%) #\N))
    (assert (= k k2))
    (call procedure trans% trans% m% n% k% alpha% a% m% b% k% z% c% m%)
    (make-matrix% m n c)))

(defun mm-hermitian% (a op-left? &optional (alpha 1))
  "Calculate alpha*A*op(A) if OP-LEFT? is NIL, and alpha*op(A)*A
  otherwise.  ALPHA defaults to 1.  op() is always conjugate
  transpose, but may be implemented as a transpose if A is real, in
  which case the two are equivalent.  This function is meant to be
  used internally, and is not exported."
  (lb-call ((common-type (lb-target-type a))
            (real-type (real-lla-type common-type))
            (procedure (lb-procedure-name common-type syrk herk))
            ((:matrix a% common-type (nrow nrow%) ncol) a)
            ((:values dim-c other-dim-a)
             (if op-left?
                 (values ncol nrow)
                 (values nrow ncol)))
            ((:output c% common-type c) (expt dim-c 2))
            ((:integer dim-c%) dim-c)
            ((:integer other-dim-a%) other-dim-a)
            ((:atom alpha% real-type) (coerce* alpha real-type))
            ((:atom beta% real-type) (zero* real-type))
            ((:char u-char%) #\U)
            ((:char op-char%) (if op-left? #\C #\N))) ; C works for real
    (call procedure u-char% op-char% dim-c% other-dim-a%
          alpha% a% nrow% beta% c% dim-c%)
    (make-matrix% dim-c dim-c c :kind :hermitian)))

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
  (declare (optimize debug))
  (bind ((diagonal-elements (elements a))
         ((:slots-read-only  nrow ncol (matrix-elements elements)) b)
         (common-type (lb-target-type diagonal-elements matrix-elements alpha))
         (alpha (coerce* alpha common-type))
         (result (make-matrix common-type nrow ncol :kind (matrix-kind b)))
         (result-elements (elements result))
         (i 0))
    (declare (fixnum i nrow ncol))
    (assert (= (length diagonal-elements) nrow) () "Dimension mismatch.")
    (with-vector-type-expansion
        (result-elements
         :other-vectors (diagonal-elements matrix-elements)
         :simple-test :other-vectors)
      (lambda (lla-type)
        `(locally 
             (declare (type ,(lla->lisp-type lla-type) alpha))
           (dotimes (col ncol)
             (dotimes (row nrow)
               (setf (aref result-elements i)
                     (* (aref matrix-elements i)
                        (aref diagonal-elements row) alpha))
               (incf i))))))
    result))

(defmethod mm ((a dense-matrix-like) (b diagonal) &optional (alpha 1))
  (declare (optimize speed))
  (bind (((:slots-read-only  nrow ncol (matrix-elements elements)) a)
         (diagonal-elements (elements b))
         (common-type (lb-target-type diagonal-elements matrix-elements alpha))
         (alpha (coerce* alpha common-type))
         (result (make-matrix common-type nrow ncol :kind (matrix-kind a)))
         (result-elements (elements result))
         (i 0))
    (assert (= (length diagonal-elements) ncol) () "Dimension mismatch.")
    (with-vector-type-expansion
        (result-elements
         :other-vectors (diagonal-elements matrix-elements)
         :simple-test :other-vectors)
      (lambda (lla-type)
        `(locally 
             (declare (type ,(lla->lisp-type lla-type) alpha))
               (dotimes (col ncol)
      (let ((d*alpha (* (aref diagonal-elements col) alpha)))
        (dotimes (row nrow)
          (setf (aref result-elements i)
                (* (aref matrix-elements i) d*alpha))
          (incf i)))))))
    result))

(defmethod mm ((a vector) (b diagonal) &optional (alpha 1))
  (mm-nv-row% a b alpha))

(defmethod mm ((a diagonal) (b vector) &optional (alpha 1))
  (mm-nv-column% a b alpha))


;;;; LU factorization

(defgeneric lu (a)
  (:documentation "LU decomposition of A"))

(defmethod lu ((a dense-matrix-like))
  (lb-call ((type (lb-target-type a))
            (procedure (lb-procedure-name type getrf))
            ((:matrix a% type (m m%) (n n%) :output lu) a)
            ((:output ipiv% :integer ipiv) (min m n))
            ((:check info%)))
    (call procedure m% n% a% m% ipiv% info%)
    (make-instance 'lu 
                   :lu-matrix (make-matrix% m n lu)
                   :ipiv ipiv)))

;;;; Hermitian factorization

(defun hermitian (a &key (component :U))
  (check-type a hermitian-matrix)
  (lb-call ((type (lb-target-type a))
            (procedure (lb-procedure-name type sytrf hetrf))
            ((:values u-char kind set-restricted?)
             (ecase component
               (:U (values #\U :upper nil))
               (:D (values #\L :lower t))))
            ((:matrix a% type (n n%) nil :output factor
                      :set-restricted? set-restricted?) a)
            ((:char u-char%) u-char)
            ((:work-queries lwork% (work% type)))
            ((:output ipiv% :integer ipiv) n)
            ((:check info%)))
    (call procedure u-char% n% a% n% ipiv% work% lwork% info%)
    (make-instance 'hermitian
                   :factor (make-matrix% n n factor :kind kind)
                   :ipiv ipiv)))

;;;; solving linear equations

(defgeneric solve (a b)
  (:documentation "Return X that solves AX=B."))

(defmethod solve (a (b vector))
  ;; Simply convert and call the method for a matrix.
  ;;
  ;; ??? check if the overhead is expensive, my hunch is that it
  ;; should be trivial -- Tamas
  (elements (solve a (as-column b))))

(defmethod solve ((a dense-matrix-like) (b dense-matrix-like))
  (lb-call ((common-type (lb-target-type a b))
            (procedure (lb-procedure-name common-type gesv))
            ((:matrix a% common-type (n n%) n2 :output :copy) a)
            ((:matrix b% common-type n3 (nrhs nrhs%) :output x) b)
            ((:work ipiv% :integer) n)  ; not saving IPIV
            ((:check info%)))
    (assert (= n n2 n3))
    (call procedure n% nrhs% a% n% ipiv% b% n% info%)
    (make-matrix% n nrhs x)))

(defmethod solve ((lu lu) (b dense-matrix-like))
  (lb-call (((:slots-r/o lu-matrix ipiv) lu)
            (common-type (lb-target-type lu-matrix b))
            (procedure (lb-procedure-name common-type getrs))
            ((:matrix lu% common-type (n n%) n2) lu-matrix)
            ((:matrix b% common-type n3 (nrhs nrhs%) :output x) b)
            ((:vector ipiv% :integer) ipiv)
            ((:char trans%) #\N)
            ((:check info%)))
    (assert (= n n2 n3))
    (call procedure trans% n% nrhs% lu% n% ipiv% b% n% info%)
    (make-matrix% n nrhs x)))

(defmethod solve ((cholesky cholesky) (b dense-matrix-like))
  (lb-call (((:slots-r/o factor) cholesky)
            (common-type (lb-target-type factor))
            (procedure (lb-procedure-name common-type potrs))
            ((:matrix factor% common-type (n n%) n2) factor)
            ((:matrix b% common-type n3 (nrhs nrhs%) :output x) b)
            ((:char u-char%) (etypecase factor
                               (upper-matrix #\U)
                               (lower-matrix #\L)))
            ((:check info%)))
    (assert (= n n2 n3))
    (call procedure u-char% n% nrhs% factor% n% b% n% info%)
    (make-matrix% n nrhs x)))

(defmethod solve ((hermitian hermitian) (b dense-matrix-like))
  (lb-call (((:slots-r/o factor ipiv) hermitian)
            (common-type (lb-target-type factor b))
            (procedure (lb-procedure-name common-type sytrs hetrs))
            ((:matrix factor% common-type (n n%) n2) factor)
            ((:matrix b% common-type n3 (nrhs nrhs%) :output x) b)
            ((:vector ipiv% :integer) ipiv)
            ((:char u-char%) (etypecase factor
                               (upper-matrix #\U)
                               (lower-matrix #\L)))
            ((:check info%)))
    (assert (= n n2 n3))
    (call procedure u-char% n% nrhs% factor% n% ipiv% b% n% info%)
    (make-matrix% n n x)))

(defun trsm% (a b side transpose-a? &optional (alpha 1))
  "Wrapper for BLAS routine xTRSM.  Calculates op(A^-1) B (if SIDE
is :LEFT) or B op(A^-1) (if SIDE is :RIGHT).  A has to be a triangular
matrix.  transpose-a? determines whether op(A) is A^T or A.  The
result is multiplied by ALPHA."
  (lb-call ((common-type (lb-target-type a b))
            (procedure (lb-procedure-name common-type trsm))
            ((:matrix a% common-type (a-n lda%) a-m) a)
            ((:matrix b% common-type (m m%) (n n%) :output result) b)
            ((:char side%) (ecase side
                             (:left (assert (= a-m m)) #\L)
                             (:right (assert (= a-n n)) #\R)))
            ((:char uplo%) (ecase (matrix-kind a)
                             (:lower #\L)
                             (:upper #\U)))
            ((:char transa%) (if transpose-a? #\C #\N))
            ((:char diag%) #\N)
            ((:atom alpha% common-type) (coerce* alpha common-type)))
    (call procedure side% uplo% transa% diag% m% n% alpha% a% lda% b% m%)
    (make-matrix% m n result)))

(defmethod solve ((a lower-matrix) (b dense-matrix-like))
  (trsm% a b :left nil))

(defmethod solve ((a upper-matrix) (b dense-matrix-like))
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
  (lb-call (((:slots-r/o lu-matrix ipiv) lu)
            (common-type (lb-target-type lu-matrix))
            (procedure (lb-procedure-name common-type getri))
            ((:matrix lu% common-type (n n%) n2 :output inverse) lu-matrix)
            ((:vector ipiv% :integer) ipiv)
            ((:work-queries lwork% (work% common-type)))
            ((:check info%)))
    (assert (= n n2))
    (call procedure n% lu% n% ipiv% work% lwork% info%)
    (make-matrix% n n inverse)))

(defmethod invert ((a hermitian-matrix) &key)
   (invert (hermitian a)))

(defmethod invert ((hf hermitian) &key)
  ;; If the FACTOR is lower triangular, we need to transpose it, as
  ;; hermitian matrices always store the upper triangle.
  (lb-call (((:slots-r/o factor ipiv) hf)
            (factor (aetypecase factor
                      (lower-matrix (transpose it))
                      (upper-matrix it)))
            (common-type (lb-target-type factor))
            (procedure (lb-procedure-name common-type sytri hetri))
            ((:matrix factor% common-type (n n%) n2 :output inverse) factor)
            ((:work work% common-type) n)
            ((:vector ipiv% :integer) ipiv)
            ((:char u-char%) #\U)
            ((:check info%)))
    (assert (= n n2))
    (call procedure u-char% n% factor% n% ipiv% work% info%)
    (make-matrix% n n inverse :kind :hermitian)))

(defun invert-triangular% (a upper? unit-diag? result-kind)
  "Invert a dense (triangular) matrix using the LAPACK routine *TRTRI.
UPPER? indicates if the matrix is in the upper or the lower triangle
of a (which needs to be a subtype of dense-matrix-like, but the type
information is not used), UNIT-DIAG? indicates whether the diagonal
is supposed to consist of 1s.  For internal use, not exported."
  (lb-call ((common-type (lb-target-type a))
            (procedure (lb-procedure-name common-type trtri))
            ((:matrix a% common-type (n n%) n2 :output inverse) a)
            ((:char u-char%) (if upper? #\U #\L))
            ((:char d-char%) (if unit-diag? #\U #\N))
            ((:check info%)))
    (assert (= n n2))
    (call procedure u-char% d-char% n% a% n% info%)
    (make-matrix% n n inverse :kind result-kind)))
  
(defmethod invert ((a upper-matrix) &key)
  (invert-triangular% a t nil :upper))

(defmethod invert ((a lower-matrix) &key)
  (invert-triangular% a nil nil :lower))

(defmethod invert ((cholesky cholesky) &key)
  ;; If the FACTOR of CHOLESKY is lower triangular, we need to
  ;; transpose it, as hermitian matrices always store the upper
  ;; triangle.
  (lb-call ((factor (aetypecase (factor cholesky)
                      (lower-matrix (transpose it))
                      (upper-matrix it)))
            (common-type (lb-target-type factor))
            (procedure (lb-procedure-name common-type potri))
            ((:matrix factor% common-type (n n%) n2 :output inverse) factor)
            ((:char u-char%) #\U)
            ((:check info%)))
    (assert (= n n2))
    (call procedure u-char% n% factor% n% info%)
    (make-matrix% n n inverse :kind :hermitian)))

(defmethod invert ((diagonal diagonal) &key (tolerance 0))
  "For pseudoinverse, suppressing diagonal elements below TOLERANCE
\(if given, otherwise / is just used without any checking."
  (make-diagonal% (emap (cond
                          ((null tolerance) #'/)
                          ((and (numberp tolerance) (<= 0 tolerance))
                           (lambda (x)
                             (if (<= (abs x) tolerance) 0 (/ x))))
                          ((and (numberp tolerance) (zerop tolerance))
                           (lambda (x)
                             (if (zerop x) 0 (/ x))))
                          (t (error "Invalid tolerance argument.")))
                        (elements diagonal))))

(defgeneric eigen (a &key vectors? check-real? &allow-other-keys)
  (:documentation "Calculate the eigenvalues and optionally the right
eigenvectors of a matrix (as columns).  Return (values eigenvalues
eigenvectors).  If check-real-p, eigenvalues of real matrices are
checked for an imaginary part and returned with the appropriate
type (compex or not).  Complex conjugate pairs of eigenvalues appear
consecutively with the eigenvalue having the positive imaginary part
first."))

(defun eigen-dense-real% (a vectors? check-real? real-type)
  "Eigenvalues and vectors for dense, real matrices."
  ;; This is probably the hairiest function in the whole library.  S/DEEV
  ;; returns real and imaginary parts for eigenvectors and eigenvalues
  ;; separately, so they have to be "zipped" together with a utility function.
  (lb-call ((procedure (lb-procedure-name real-type geev))
            ((:matrix a% real-type (n n%) n2 :output :copy) a)
            ((:output w% real-type w) (* 2 n))
            ((:output vr% real-type vr) (if vectors? (expt n 2) 0))
            ((:work-queries lwork% (work% real-type)))
            ((:char n-char%) #\N)
            ((:char v-char%) #\V)
            ((:check info%)))
    (assert (= n n2))
    (call procedure n-char% v-char% n% a% n% 
          w% (inc-pointer w% (* n (foreign-size* real-type))) ; eigenvalues
          (null-pointer) n% vr% n%                            ; eigenvectors
          work% lwork% info%)
    (zip-eigen% w (when vectors? vr) check-real? real-type)))

(defun eigen-dense-complex% (a vectors? complex-type)
  "Eigenvalues and vectors for dense, complex matrices."
  (lb-call ((procedure (lb-procedure-name complex-type geev))
            ((:matrix a% complex-type (n n%) n2 :output :copy) a)
            ((:output w% complex-type w) n)
            ((:work rwork% complex-type) n)
            ((:char n-char%) #\N)
            ((:char v-char%) #\V)
            ((:output vr% complex-type vr) (if vectors? (expt n 2) 0))
            ((:work-queries lwork% (work% complex-type)))
            ((:check info%)))
    (assert (= n n2))
    (call procedure n-char% v-char% n% a% n% w%
          (null-pointer) n% vr% n% work% lwork% rwork% info%)
    (values w (when vectors?
                (make-matrix% n n vr)))))

(defmethod eigen ((a dense-matrix-like) &key vectors? check-real?)
  ;; Unfortunately, the LAPACK interface for real and complex cases is
  ;; different.
  (let ((type (lb-target-type a)))
    (ecase type
      ((:single :double)
         (eigen-dense-real% a vectors? check-real? type))
      ((:complex-single :complex-double)
         (eigen-dense-complex% a vectors? type)))))


(defmethod eigen ((a hermitian-matrix) &key vectors? check-real?)
  ;; Uses simple driver, where "simple" means "silly collection of
  ;; special cases for all combinations".
  (declare (ignore check-real?))        ; eigenvalues are always real
  (lb-call ((common-type (lb-target-type a))
            (complex? (lla-complex? common-type))
            (real-type (real-lla-type common-type))
            (procedure (lb-procedure-name common-type syev heev))
            ((:char nv-char%) (if vectors? #\V #\N))
            ((:char u-char%) #\U)       ; upper triangle
            ((:matrix a% common-type (n n%) n2 :output v :set-restricted? nil) a)
            ((:output w% real-type w) n) ; eigenvalues
            ((:work rwork% real-type) (if complex? (max 1 (- (* 3 n) 2)) 0))
            ((:work-queries lwork% (work% real-type)))
            ((:check info%)))
    (assert (= n n2))
    (if complex?
        (call procedure nv-char% u-char% n% a% n% w%
              work% lwork% rwork% info%)
        (call procedure nv-char% u-char% n% a% n% w%
              work% lwork% info%))
    (values w (when vectors? (make-matrix% n n v)))))

;;;;
;;;; least squares calculations
;;;;

(defgeneric least-squares (y x &key method &allow-other-keys)
  (:documentation "Return (values b ss nu other-values), where beta = argmin_b
L2norm( y-Xb ), solving a least squares problem, SS is the sum of squares for
each column of Y, and NU is the degrees of freedom.  Y can have multiple
columns, in which case X will have the same number of columns, each
corresponding to a different column of Y.  METHOD selects the method to call,
methods may accept other keyword arguments and return additional values as a
list property list in OTHER-VALUES."))

(define-condition not-enough-columns (error)
  ())

(defmethod least-squares ((y dense-matrix-like) (x dense-matrix-like) &rest
                          named-pairs &key (method :svd-d) &allow-other-keys)
  (apply (ecase method
           (:qr #'least-squares-qr)
           (:svd-d #'least-squares-svd-d))
         y x 
         :allow-other-keys t
         named-pairs))

;; ;;; univariate versions of least squares: vector ~ vector, vector ~ matrix

(defmethod least-squares ((y vector) (x dense-matrix-like) &rest named-pairs)
  (bind (((:values b ss nu other-values)
          (apply #'least-squares (as-column y) x named-pairs)))
    (values (elements b) (aref ss 0) nu other-values)))

(defmethod least-squares ((y vector) (x vector) &rest named-pairs)
  (bind (((:values b ss nu other-values)
          (apply #'least-squares (as-column y) (as-column x) named-pairs)))
    (values (aref (elements b) 0) (aref ss 0) nu other-values)))

(defun least-squares-qr (y x &key &allow-other-keys)
  "Least squares using QR decomposition.  Additional values returned: the QR
decomposition of X.  See LEAST-SQUARES for additional documentation.  Usage
note: SVD-based methods are recommended over this one, unless X is
well-conditioned."
  ;; Note: the naming convention (y,X,b) is different from LAPACK's
  ;; (b,A,x).  Sorry if this creates confusion, I decided to follow
  ;; standard statistical notation.
  (lb-call ((common-type (lb-target-type x y))
            (real-type (real-lla-type common-type))
            (procedure (lb-procedure-name common-type gels))
            ((:matrix x% common-type (m m%) (n n%) :output qr) x)
            ((:matrix y% common-type m2 (nrhs nrhs%) :output b) y)
            ((:char n-char%) #\N)
            ((:work-queries lwork% (work% real-type)))
            ((:check info%)))
    (assert (= m m2))
    (unless (<= n m)
      (error 'not-enough-columns))
    (call procedure n-char% m% n% nrhs% x% m% y% m% work% lwork% info%)
    (values 
      (matrix-from-first-rows b n nrhs m)
      (sum-last-rows b m nrhs n)
      (- m n)
      `(:qr ,(make-instance 'qr :qr-matrix (make-matrix% m n qr))))))

(defun least-squares-svd-d (y x &key (rcond -1))
  (bind ((common-type (lb-target-type x y))
         (real-type (real-lla-type common-type))
         (procedure (lb-procedure-name common-type gelsd)))
    (lb-call (((:matrix x% common-type (m m%) (n n%) :output :copy) x)
              ((:matrix y% common-type m2 (nrhs nrhs%) :output b) y)
              ((:work-queries lwork% 
                              (work% common-type)
                              (rwork% real-type t)
                              (iwork% :integer)))
              ((:check info%))
              ((:output s% real-type s) (min m n))
              ((:output rank% :integer rank) 1)
              ((:atom rcond% real-type) (coerce* rcond real-type)))
      (assert (= m m2))
      (unless (<= n m)
        (error 'not-enough-columns))
      (if (lla-complex? common-type)
          (call procedure m% n% nrhs% x% m% y% m% s% rcond% rank% work% lwork% 
                rwork% iwork% info%)
          (call procedure m% n% nrhs% x% m% y% m% s% rcond% rank% work% lwork% 
                iwork% info%))
      (values 
        (matrix-from-first-rows b n nrhs m)
        (sum-last-rows b m nrhs n)
        (- m n)
        `(:s ,s :rank ,(aref rank 0))))))

(defun qr-xx-inverse (qr)
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

(defmethod constrained-least-squares ((y vector) (x dense-matrix-like)
                                      (z dense-matrix-like) (w vector))
  "Solve the (linearly) constrained least squares problem min_b |y-Xb|_2
  subject to Zx=w."
  ;; Note: mapping between the function parameters/variables and
  ;; LAPACK counterparts is as follows: y->c, X->A, x->b, Z->B, w->d,
  (lb-call ((common-type (lb-target-type y x z w))
            (procedure (lb-procedure-name common-type gglse))
            ((:matrix x% common-type (m m%) (n n%) :output :copy) x)
            ((:matrix z% common-type (p p%) n2 :output :copy) z)
            ((:vector y% common-type :copy) y)
            ((:vector w% common-type :copy) w)
            ((:output b% common-type b) n)
            ((:work-queries lwork% (work% common-type)))
            ((:check info%)))
    ;; !! after call, x and z contain decompositions, currently not collected
    ;; !! after call, y contains sum of squares, currently not collected
    (assert (= n n2) () "Dimension mismatch between z and x")
    (assert (= m (length y)) () "Dimension mismatch between x and y")
    (assert (= p (length w)) () "Dimension mismatch between z and w")
    (call procedure m% n% p% x% m% z% p% y% w% b% work% lwork% info%)
    b))


;; ;;;; Cholesky factorization

(defgeneric cholesky (a &optional component)
  (:documentation "Cholesky factorization.  Component is :L or :U (default)."))

(defmethod cholesky ((a hermitian-matrix) &optional (component :U))
  (lb-call ((common-type (lb-target-type a))
            (procedure (lb-procedure-name common-type potrf))
            ((:values u-char kind) (ecase component
                                     (:U (values #\U :upper))
                                     (:L (values #\L :lower))))
            ((:matrix a% common-type (n n%) n2 :output cholesky) a)
            ((:char u-char%) u-char)
            ((:check info%)))
    (assert (= n n2))
    (call procedure u-char% n% a% n% info%)
    (make-instance 'cholesky 
                   :factor (make-matrix% n n cholesky :kind kind))))

(defmethod reconstruct ((mf cholesky))
  (bind (((:slots-read-only factor) mf))
    (etypecase factor
      (upper-matrix (mm t factor))
      (lower-matrix (mm factor t)))))


;; ;;; SVD

(defgeneric svd (a &key left right)
  (:documentation "Return singular value decomposition A, with left- and right
singular vectors as second and third values (when requested, see vector
specifications).  Valid vector specifications are :NONE, :SINGULAR (singular
vectors only) and :ALL.  Return values S (singular values, descending order, as
a DIAGONAL), U (left singular vectors, DENSE-MATRIX, or NIL), VT ([conjugate]
transpose of right singular vectors, DENSE-MATRIX, or NIL)."))

(defmethod svd ((a dense-matrix-like) &key (left :none) (right :none))
  (lb-call ((type (lb-target-type a))
            (real-type (real-lla-type type))
            (complex? (lla-complex? type))
            (procedure (lb-procedure-name type gesvd))
            ((:matrix a% type (m m%) (n n%) :output :copy) a)
            (min-mn (min m n))
            ((:values u-ncol jobu) (ecase left
                                     (:none (values 1 #\N)) ; LAPACK needs >=1
                                     (:singular (values min-mn #\S))
                                     (:all (values m #\A))))
            ((:values vt-nrow jobvt) (ecase right
                                       (:none (values 1 #\N)) ; LAPACK needs >=1
                                       (:singular (values min-mn #\S))
                                       (:all (values n #\A))))
            ((:char jobu%) jobu)
            ((:char jobvt%) jobvt)
            ((:integer vt-nrow%) vt-nrow)
            ((:output s% real-type s) min-mn)
            ((:output u% type u) (* m u-ncol))
            ((:output vt% type vt) (* vt-nrow n))
            ((:work rwork% real-type) (if complex? (* 5 min-mn) 0))
            ((:work-queries lwork% (work% type)))
            ((:check info%)))
    (with-lapack-traps-masked
      (if complex?
          (call procedure jobu% jobvt% m% n% a% m% s% u% m% vt% vt-nrow%
                work% lwork% rwork% info%)
          (call procedure jobu% jobvt% m% n% a% m% s% u% m% vt% vt-nrow%
                work% lwork% info%)))
    (values (make-diagonal% s)
            (unless (eq left :none) (make-matrix% m u-ncol u))
            (unless (eq right :none) (make-matrix% vt-nrow n vt)))))

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
         ((:values matrix squared?)
          (cond 
            ((aand logrc-threshold (< ratio (- it))) ; fewer rows than columns
             (values (mm matrix t) t))
            ((aand logrc-threshold (< it ratio)) ; more rows than columns
             (values (mm t matrix) t))
            (t                          ; no matrix multiplication
             (values matrix nil))))
         (d (svd matrix))
         (d-elements (elements d)))
    (map-into d-elements (if squared? (lambda (x) (sqrt (abs x))) #'abs)
              d-elements)
    (let ((threshold (aif threshold
                          it
                          (* (max nrow ncol) (epsilon* (lb-target-type d))
                             (reduce #'max d-elements)))))
      (values (count-if (lambda (x) (<= threshold (abs x))) (elements d)) d))))

;;; determinants

(defgeneric logdet (matrix)
  (:documentation "Logarithm of the determinant of a matrix."))

(defun det (matrix)
  "Determinant of a matrix.  If you need the log of this, use LOGDET
  directly."
  (exp (logdet matrix)))

(defun diagonal-log-sum% (matrix)
  "Sum of the log of the elements in the diagonal."
  (assert (square-matrix? matrix))
  (bind (((:lla-matrix matrix :nrow nrow) matrix))
    (iter
      (for i :from 0 :below nrow)
      (summing (log (matrix (matrix-index i i)))))))

(defmethod logdet ((matrix lower-matrix))
  (diagonal-log-sum% matrix))

(defmethod logdet ((matrix upper-matrix))
  (diagonal-log-sum% matrix))

(defmethod logdet ((matrix hermitian-matrix))
  (reduce #'+ (eigen matrix) :key #'log))

(defun matrix-cond (matrix)
  "Calculate the condition number of a matrix (with respect to the Euclidean
  norm).  Implemented as the ratio of the largest and smallest singular values."
  (check-type matrix dense-matrix-like)
  (let* ((s (elements (svd matrix)))
         (s-min (aref s (1- (length s)))))
    (if (zerop s-min)
        (error "Matrix is not full rank.")
        (/ (aref s 0) s-min))))
