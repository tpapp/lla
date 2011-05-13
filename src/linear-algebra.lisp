(in-package :lla)

;;; Higher level linear algebra functions defined here.
;;;
;;; General convention for vectors in places of matrices: should be
;;; interpreted as a conforming vector.  If the result is an 1xn or nx1
;;; matrix, it should be converted to a vector iff some related argument was a
;;; vector.  For example, (solve a b) should be a vector iff b is a vector,
;;; otherwise it should be a matrix.  In case of ambiguity (eg (mm a b) could
;;; be both (dot a b) or (outer a b) when a and b are vectors) signal an
;;; error.

;;; General notes about LAPACKE/CBLAS:
;;; 
;;; The leading dimension for row-major matrices is the SECOND dimension
;;; (number of columns).

;;; Helper functions

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

(defun mm-internal (a b &key a-op b-op (alpha 1))
  (lb-call (((:values a0 a1 a-orientation lda a-op)
             (maybe-vector-as-matrix a :row a-op))
            ((:values b0 b1 b-orientation ldb b-op)
             (maybe-vector-as-matrix b :column b-op))
            (common-type (common-float-type a b alpha))
            ((:blas gemm common-type))
            ((:array a% common-type) a)
            ((:array b% common-type) b)
            ((:output c% common-type c)
             (cond
               ((and a-orientation b-orientation) 
                (error "Can't use MM on two vectors, see DOT or OUTER."))
               (a-orientation b1)
               (b-orientation a0)
               (t (list a0 b1)))))
    (assert (= a1 b0))
    (call a-op b-op a0 b1 a1 (coerce* alpha common-type) a% lda b% ldb
          (zero* common-type) c% b1)
    c))

(defgeneric mm (a b &optional alpha)
  (:documentation
   "multiply A and B, also by the scalar alpha (defaults to 1)."))

(defun mmm (&rest matrices)
  (reduce #'mm matrices))

(defmethod mm ((a array) (b array) &optional (alpha 1))
  (mm-internal a b :alpha alpha))

;;; !! this is how we could speed things up with compiler macros: have a
;;; !! compiler macro transform (MM (TRANSPOSE FOO) BAR) to (MM-TN FOO BAR),
;;; !! which is a generic function and would use MM-INTERNAL accordingly.  Or,
;;; !! we could have (lazy-transpose* FOO) explicitly, and have MM dispatch on
;;; !! that.  Will investigate.

(defun mm-hermitian% (a op-left? &optional (alpha 1))
  "Calculate alpha*A*op(A) if OP-LEFT? is NIL, and alpha*op(A)*A
  otherwise.  ALPHA defaults to 1.  op() is always conjugate
  transpose, but may be implemented as a transpose if A is real, in
  which case the two are equivalent.  This function is meant to be
  used internally, and is not exported."
  (lb-call ((common-type (common-float-type a))
            (real-type (real-lla-type common-type))
            ((:blas (syrk herk) common-type))
            ((:array a% common-type :dimensions (nrow ncol)) a)
            ((:values dim-c other-dim-a)
             (if op-left?
                 (values ncol nrow)
                 (values nrow ncol)))
            ((:output c% common-type c) (list dim-c dim-c)))
    (call (hermitian-orientation :blas) (lb-transpose op-left? :blas t)
          dim-c other-dim-a (coerce* alpha real-type) a% ncol
          (zero* real-type) c% dim-c)
    (make-instance 'hermitian-matrix :elements c)))

(defmethod mm ((a array) (b (eql t)) &optional (alpha 1))
  ;; A A^T
  (mm-hermitian% a nil alpha))

(defmethod mm ((a (eql t)) (b array) &optional (alpha 1))
  ;; A^T A
  (mm-hermitian% b t alpha))

(defmethod mm ((a wrapped-matrix) (b wrapped-matrix) &optional (alpha 1))
  (mm (as-array a) (as-array b) alpha))

(defmethod mm ((a wrapped-matrix) b &optional (alpha 1))
  (mm (as-array a) b alpha))

(defmethod mm (a (b wrapped-matrix) &optional (alpha 1))
  (mm a (as-array b) alpha))

;; ;;; mm for diagonal matrices

;; (defmethod mm ((a diagonal) (b dense-matrix-like) &optional (alpha 1))
;;   (declare (optimize speed))
;;   (bind ((diagonal-elements (elements a))
;;          ((:slots-read-only  nrow ncol (matrix-elements elements)) b)
;;          (common-type (lb-target-type diagonal-elements matrix-elements alpha))
;;          (alpha (coerce* alpha common-type))
;;          (result (make-matrix nrow ncol common-type :kind (matrix-kind b)))
;;          (result-elements (elements result))
;;          (i 0))
;;     (declare (fixnum i nrow ncol))
;;     (assert (= (length diagonal-elements) nrow) () "Dimension mismatch.")
;;     (with-vector-type-expansion
;;         (result-elements
;;          :other-vectors (diagonal-elements matrix-elements)
;;          :simple-test :other-vectors)
;;       (lambda (lla-type)
;;         `(locally 
;;              (declare (type ,(lla->lisp-type lla-type) alpha))
;;            (dotimes (col ncol)
;;              (dotimes (row nrow)
;;                (setf (aref result-elements i)
;;                      (* (aref matrix-elements i)
;;                         (aref diagonal-elements row) alpha))
;;                (incf i))))))
;;     result))

;; (defmethod mm ((a dense-matrix-like) (b diagonal) &optional (alpha 1))
;;   (declare (optimize speed))
;;   (bind (((:slots-read-only  nrow ncol (matrix-elements elements)) a)
;;          (diagonal-elements (elements b))
;;          (common-type (lb-target-type diagonal-elements matrix-elements alpha))
;;          (alpha (coerce* alpha common-type))
;;          (result (make-matrix nrow ncol common-type :kind (matrix-kind a)))
;;          (result-elements (elements result))
;;          (i 0))
;;     (assert (= (length diagonal-elements) ncol) () "Dimension mismatch.")
;;     (with-vector-type-expansion
;;         (result-elements
;;          :other-vectors (diagonal-elements matrix-elements)
;;          :simple-test :other-vectors)
;;       (lambda (lla-type)
;;         `(locally 
;;              (declare (type ,(lla->lisp-type lla-type) alpha))
;;                (dotimes (col ncol)
;;       (let ((d*alpha (* (aref diagonal-elements col) alpha)))
;;         (dotimes (row nrow)
;;           (setf (aref result-elements i)
;;                 (* (aref matrix-elements i) d*alpha))
;;           (incf i)))))))
;;     result))

;; (defmethod mm ((a vector) (b diagonal) &optional (alpha 1))
;;   (mm-row% a b alpha))

;; (defmethod mm ((a diagonal) (b vector) &optional (alpha 1))
;;   (mm-column% a b alpha))

;;; hermitian (symmetric) updates

;; (defun update-hermitian (a x &optional (alpha 1))
;;   "Return alpha x x^H + A.  A has to be a hermitian matrix, and ALPHA real."
;;   (check-type a hermitian-matrix)
;;   (check-type x vector)
;;   (lb-call ((common-type (lb-target-type a x alpha))
;;             (real-type (real-lla-type common-type))
;;             (procedure (lb-procedure-name common-type syr her))
;;             ((:matrix a% common-type (n n%) n2 :output c) a)
;;             ((:vector x% common-type n3) x)
;;             ((:atom alpha% real-type) (coerce* alpha real-type))
;;             ((:char u%) #\U)
;;             ((:integer one%) 1))
;;     (assert (= n n2 n3))
;;     (call procedure u% n% alpha% x% one% a% n%)
;;     (make-matrix% n n c :kind :hermitian)))

;; (defun update-hermitian2 (a x y &optional (alpha 1))
;;   "Return alpha x y^H + y (alpha x)^H + A.  A has to be a hermitian matrix."
;;   (check-type a hermitian-matrix)
;;   (check-type x vector)
;;   (check-type y vector)
;;   (lb-call ((common-type (lb-target-type a x y alpha))
;;             (procedure (lb-procedure-name common-type syr2 her2))
;;             ((:matrix a% common-type (n n%) n2 :output c) a)
;;             ((:vector x% common-type n3) x)
;;             ((:vector y% common-type n4) y)
;;             ((:atom alpha% common-type) (coerce* alpha common-type))
;;             ((:char u%) #\U)
;;             ((:integer one%) 1))
;;     (assert (= n n2 n3 n4))
;;     (call procedure u% n% alpha% x% one% y% one% a% n%)
;;     (make-matrix% n n c :kind :hermitian)))


;;;; LU factorization

(defgeneric lu (a)
  (:documentation "LU decomposition of A"))

(defmethod lu ((a array))
  (lb-call ((type (common-float-type a))
            ((:lapack getrf type))
            ((:array a% type :dimensions (m n) :output lu) a)
            ((:output ipiv% :integer ipiv) (min m n)))
    (call m n a% n ipiv%)
    (make-instance 'lu :lu lu :ipiv ipiv)))

;;;; Hermitian factorization

(defgeneric hermitian-factorization (a)
  (:documentation "Compute the hermitian factorization."))

(defmethod hermitian-factorization ((a hermitian-matrix))
  (lb-call ((a (elements a))
            (type (common-float-type a))
            ((:lapack (sytrf hetrf) type))
            ((:array a% type :dimensions (n n2) :output factor) a)
            ((:output ipiv% :integer ipiv) n))
    (assert (= n n2))
    (call (hermitian-orientation :lapack) n a% n ipiv%)
    (make-instance 'hermitian-factorization :factor factor :ipiv ipiv)))

;;;; solving linear equations

(defgeneric solve (a b)
  (:documentation "Return X that solves AX=B.  When B is a vector, so is X."))

(defmethod solve (a (b wrapped-matrix))
  (solve a (as-array b)))

(defmethod solve ((lu lu) (b array))
  (lb-call (((:values b0 b1 b-orientation ldb)
             (maybe-vector-as-matrix b :column))
            ((:slots-r/o lu ipiv) lu)
            (common-type (common-float-type lu b))
            ((:lapack getrs common-type))
            ((:array lu% common-type :dimensions (lu0 lu1)) lu)
            ((:array b% common-type :output x :output-dimensions 
                     (vector-or-matrix-dimensions lu1 b1 b-orientation)) b)
            ((:array ipiv% :integer :dimensions (i0)) ipiv))
    (assert (= lu0 lu1 b0 i0))
    (call +n+ lu0 b1 lu% lu1 ipiv% b% ldb)
    x))

(defmethod solve ((a array) (b array))
  (lb-call (((:values b0 b1 b-orientation ldb)
             (maybe-vector-as-matrix b :column))
            (common-type (common-float-type a b))
            ((:lapack gesv common-type))
            ((:array a% common-type :dimensions (a0 a1)) a)
            ((:array b% common-type :output x :output-dimensions
                     (vector-or-matrix-dimensions a1 b1 b-orientation)) b)
            ((:work ipiv% :integer) a0)) ; not saving ipiv
    (assert (= a0 a1 b0))
    (call a0 b1 a% a1 ipiv% b% ldb)
    x))

;; (defmethod solve ((cholesky cholesky) (b dense-matrix-like))
;;   (lb-call (((:slots-r/o factor) cholesky)
;;             (common-type (common-float-type factor))
;;             (procedure (lb-procedure-name common-type potrs))
;;             ((:matrix factor% common-type (n n%) n2) factor)
;;             ((:matrix b% common-type n3 (nrhs nrhs%) :output x) b)
;;             ((:char u-char%) (etypecase factor
;;                                (upper-matrix #\U)
;;                                (lower-matrix #\L)))
;;             ((:check info%)))
;;     (assert (= n n2 n3))
;;     (call procedure u-char% n% nrhs% factor% n% b% n% info%)
;;     (make-matrix% n nrhs x)))

(defmethod solve ((hermitian-factorization hermitian-factorization) (b array))
  (lb-call (((:slots-r/o factor ipiv) hermitian-factorization)
            (common-type (common-float-type factor b))
            ((:lapack (sytrs hetrs) common-type))
            ((:array factor% common-type :dimensions (n n2)) factor)
            ((:values n3 nrhs b-orientation ldb)
             (maybe-vector-as-matrix b :column))
            ((:array b% common-type :output x :output-dimensions
                      (vector-or-matrix-dimensions n nrhs b-orientation)) b)
            ((:array ipiv% :integer :dimensions n4) ipiv))
    (assert (= n n2 n3 n4))
    (call (hermitian-orientation :lapack) n nrhs factor% n ipiv% b% ldb)
    x))

(defmethod solve ((hermitian-matrix hermitian-matrix) (b array))
  (solve (hermitian-factorization hermitian-matrix) b))

(defun trsm% (a a-uplo b side transpose-a? &optional (alpha 1))
  "Wrapper for BLAS routine xTRSM.  Calculates op(A^-1) B (if SIDE
is :LEFT) or B op(A^-1) (if SIDE is :RIGHT).  A has to be a triangular
matrix.  transpose-a? determines whether op(A) is A^T or A.  The
result is multiplied by ALPHA."
  (lb-call ((common-type (common-float-type a b))
            ((:blas trsm common-type))
            ((:values b0 b1 b-orientation ldb)
             (maybe-vector-as-matrix b :column))
            ((:array a% common-type :dimensions (a0 a1)) a)
            ((:array b% common-type :output result :output-dimensions
                     (vector-or-matrix-dimensions b0 b1 b-orientation)) b)
            (side (ecase side
                    (:left (assert (= a0 b0)) :CBLASLEFT)
                    (:right (assert (= a0 b1)) :CBLASRIGHT)))
            (trans (lb-transpose transpose-a? :blas t))
            (alpha (coerce* alpha common-type)))
    (assert (= a0 a1))
    (call side a-uplo trans :CBLASNONUNIT b0 b1 alpha a% a0 b% ldb)
    result))

(defmethod solve ((a lower-triangular-matrix) b)
  (trsm% (elements a) :CBLASLOWER (as-array b) :left nil))

(defmethod solve ((a upper-triangular-matrix) b)
  (trsm% (elements a) :CBLASUPPER (as-array b) :left nil))

;; (defmethod solve ((a diagonal) (b dense-matrix-like))
;;   (mm (e/ a) b))

(defgeneric invert (a &key &allow-other-keys)
  (:documentation "Invert A.  The inverse of matrix factorizations are other
  factorizations when appropriate, otherwise the result is a matrix.  Usage
  note: inverting dense matrices is unnecessary and unwise in most cases,
  because it is numerically unstable.  If you are solving many Ax=b equations
  with the same A, use a matrix factorization like LU."))

(defmethod invert ((a array) &key) (invert (lu a)))

(defmethod invert ((lu lu) &key)
  (lb-call (((:slots-r/o lu ipiv) lu)
            (common-type (common-float-type lu))
            ((:lapack getri common-type))
            ((:array lu% common-type :dimensions (n n2) :output inverse) lu)
            ((:array ipiv% :integer :dimensions n3) ipiv))
    (assert (= n n2 n3))
    (call n lu% n ipiv%)
    inverse))

(defmethod invert ((hermitian-factorization hermitian-factorization) &key)
  (lb-call (((:slots-r/o factor ipiv) hermitian-factorization)
            (common-type (common-float-type factor))
            ((:lapack (sytri hetri) common-type))
            ((:array factor% common-type :dimensions (n n2) :output inverse)
             factor)
            ((:array ipiv% :integer :dimensions n3) ipiv))
    (assert (= n n2 n3))
    (call (hermitian-orientation :lapack) n factor% n ipiv%)
    (make-instance 'hermitian-matrix :elements inverse)))

(defmethod invert ((a hermitian-matrix) &key)
  (invert (hermitian-factorization a)))

(defun invert-triangular% (a upper? unit-diag? kind)
  "Invert a dense (triangular) matrix using the LAPACK routine *TRTRI.
UPPER? indicates if the matrix is in the upper or the lower triangle of a
matrix, UNIT-DIAG? indicates whether the diagonal is supposed to consist of
1s.  For internal use, not exported."
  (lb-call ((common-type (common-float-type a))
            ((:lapack trtri common-type))
            ((:array a% common-type :dimensions (n n2) :output inverse) a)
            (u-char (if upper? +u+ +l+))
            (d-char (if unit-diag? +u+ +n+)))
    (assert (= n n2))
    (call u-char d-char n a% n)
    (make-matrix kind (list n n) :initial-contents inverse)))
  
(defmethod invert ((a upper-triangular-matrix) &key)
  (invert-triangular% (elements a) t nil :upper))

(defmethod invert ((a lower-triangular-matrix) &key)
  (invert-triangular% (elements a) nil nil :lower))

;; (defmethod invert ((cholesky cholesky) &key)
;;   ;; If the FACTOR of CHOLESKY is lower triangular, we need to
;;   ;; transpose it, as hermitian matrices always store the upper
;;   ;; triangle.
;;   (lb-call ((factor (aetypecase (factor cholesky)
;;                       (lower-matrix (transpose it))
;;                       (upper-matrix it)))
;;             (common-type (common-float-type factor))
;;             (procedure (lb-procedure-name common-type potri))
;;             ((:matrix factor% common-type (n n%) n2 :output inverse) factor)
;;             ((:char u-char%) #\U)
;;             ((:check info%)))
;;     (assert (= n n2))
;;     (call procedure u-char% n% factor% n% info%)
;;     (make-matrix% n n inverse :kind :hermitian)))

;; (defmethod invert ((diagonal diagonal) &key (tolerance 0))
;;   "For pseudoinverse, suppressing diagonal elements below TOLERANCE
;; \(if given, otherwise / is just used without any checking."
;;   (make-diagonal% (emap (cond
;;                           ((null tolerance) #'/)
;;                           ((and (numberp tolerance) (<= 0 tolerance))
;;                            (lambda (x)
;;                              (if (<= (abs x) tolerance) 0 (/ x))))
;;                           ((and (numberp tolerance) (zerop tolerance))
;;                            (lambda (x)
;;                              (if (zerop x) 0 (/ x))))
;;                           (t (error "Invalid tolerance argument.")))
;;                         (elements diagonal))))

;; (defgeneric eigen (a &key vectors? check-real? &allow-other-keys)
;;   (:documentation "Calculate the eigenvalues and optionally the right
;; eigenvectors of a matrix (as columns).  Return (values eigenvalues
;; eigenvectors).  If check-real-p, eigenvalues of real matrices are
;; checked for an imaginary part and returned with the appropriate
;; type (compex or not).  Complex conjugate pairs of eigenvalues appear
;; consecutively with the eigenvalue having the positive imaginary part
;; first."))

;; (defun eigen-dense-real% (a vectors? check-real? real-type)
;;   "Eigenvalues and vectors for dense, real matrices."
;;   ;; This is probably the hairiest function in the whole library.  S/DEEV
;;   ;; returns real and imaginary parts for eigenvectors and eigenvalues
;;   ;; separately, so they have to be "zipped" together with a utility function.
;;   (lb-call ((procedure (lb-procedure-name real-type geev))
;;             ((:matrix a% real-type (n n%) n2 :output :copy) a)
;;             ((:output w% real-type w) (* 2 n))
;;             ((:output vr% real-type vr) (if vectors? (expt n 2) 0))
;;             ((:work-queries lwork% (work% real-type)))
;;             ((:char n-char%) #\N)
;;             ((:char v-char%) #\V)
;;             ((:check info%)))
;;     (assert (= n n2))
;;     (call procedure n-char% v-char% n% a% n% 
;;           w% (inc-pointer w% (* n (foreign-size* real-type))) ; eigenvalues
;;           (null-pointer) n% vr% n%                            ; eigenvectors
;;           work% lwork% info%)
;;     (zip-eigen% w (when vectors? vr) check-real? real-type)))

;; (defun eigen-dense-complex% (a vectors? complex-type)
;;   "Eigenvalues and vectors for dense, complex matrices."
;;   (lb-call ((procedure (lb-procedure-name complex-type geev))
;;             ((:matrix a% complex-type (n n%) n2 :output :copy) a)
;;             ((:output w% complex-type w) n)
;;             ((:work rwork% complex-type) n)
;;             ((:char n-char%) #\N)
;;             ((:char v-char%) #\V)
;;             ((:output vr% complex-type vr) (if vectors? (expt n 2) 0))
;;             ((:work-queries lwork% (work% complex-type)))
;;             ((:check info%)))
;;     (assert (= n n2))
;;     (call procedure n-char% v-char% n% a% n% w%
;;           (null-pointer) n% vr% n% work% lwork% rwork% info%)
;;     (values w (when vectors?
;;                 (make-matrix% n n vr)))))

;; (defmethod eigen ((a dense-matrix-like) &key vectors? check-real?)
;;   ;; Unfortunately, the LAPACK interface for real and complex cases is
;;   ;; different.
;;   (let ((type (common-float-type a)))
;;     (ecase type
;;       ((:single :double)
;;          (eigen-dense-real% a vectors? check-real? type))
;;       ((:complex-single :complex-double)
;;          (eigen-dense-complex% a vectors? type)))))


;; (defmethod eigen ((a hermitian-matrix) &key vectors? check-real?)
;;   ;; Uses simple driver, where "simple" means "silly collection of
;;   ;; special cases for all combinations".
;;   (declare (ignore check-real?))        ; eigenvalues are always real
;;   (lb-call ((common-type (common-float-type a))
;;             (complex? (lla-complex? common-type))
;;             (real-type (real-lla-type common-type))
;;             (procedure (lb-procedure-name common-type syev heev))
;;             ((:char nv-char%) (if vectors? #\V #\N))
;;             ((:char u-char%) #\U)       ; upper triangle
;;             ((:matrix a% common-type (n n%) n2 :output v :set-restricted? nil) a)
;;             ((:output w% real-type w) n) ; eigenvalues
;;             ((:work rwork% real-type) (if complex? (max 1 (- (* 3 n) 2)) 0))
;;             ((:work-queries lwork% (work% real-type)))
;;             ((:check info%)))
;;     (assert (= n n2))
;;     (if complex?
;;         (call procedure nv-char% u-char% n% a% n% w%
;;               work% lwork% rwork% info%)
;;         (call procedure nv-char% u-char% n% a% n% w%
;;               work% lwork% info%))
;;     (values w (when vectors? (make-matrix% n n v)))))

;;;;
;;;; least squares calculations
;;;;
;;;; All least squares functions return (values b ss nu ...), where beta =
;;;; argmin_b L2norm( y-Xb ), solving a least squares problem, SS is the sum of
;;;; squares for each column of Y, and NU is the degrees of freedom.  Y can have
;;;; multiple columns, in which case X will have the same number of columns, each
;;;; corresponding to a different column of Y.

(defun least-squares (y x &rest rest &key (method :qr) &allow-other-keys)
  (ecase method
    (:qr (apply #'least-squares-qr y x rest))))

(define-condition not-enough-columns (error)
  ())

(defun last-rows-ss (matrix nrhs common-type)
  "Calculate the sum of squares of the last rows of MATRIX columnwise,
omitting the first NRHS rows.  Used for interfacing with xGELS."
  (bind (((m n) (array-dimensions matrix))
         (real-type (real-lla-type common-type))
         (sum (make-array* n real-type (zero* real-type)))
         (matrix-index (array-row-major-index matrix nrhs 0)))
    (loop repeat (- m nrhs) do
      (dotimes (sum-index n)
        (incf (aref sum sum-index)
              (conjugate-square (row-major-aref matrix matrix-index)))
        (incf matrix-index)))
    sum))

(defun least-squares-qr (y x &key &allow-other-keys)
  "Least squares using QR decomposition.  Additional values returned: the QR
decomposition of X.  See LEAST-SQUARES for additional documentation.  Usage
note: SVD-based methods are recommended over this one, unless X is
well-conditioned."
  ;; Note: the naming convention (y,X,b) is different from LAPACK's
  ;; (b,A,x).  Sorry if this creates confusion, I decided to follow
  ;; standard statistical notation.
  (lb-call ((common-type (common-float-type x y))
            ((:values m n x-orientation ldx)
             (maybe-vector-as-matrix x :column))
            ((:values m2 nrhs y-orientation ldy)
             (maybe-vector-as-matrix y :column))
            ((:lapack gels common-type))
            ((:array x% common-type
                     :output qr :output-dimensions (list m n)) x)
            ((:array y% common-type
                     :output b :output-dimensions (list m nrhs)) y))
    (assert (= m m2))
    (unless (<= n m)
      (error 'not-enough-columns))
    (call +n+ m n nrhs x% ldx y% ldy)
    (values
      (maybe-pick-first-element (matrix-from-first-rows b n y-orientation)
                                x-orientation)
      (maybe-pick-first-element (last-rows-ss b n common-type)
                                y-orientation)
      (- m n)
      (make-instance 'qr :qr qr))))

;; (defun least-squares-svd-d (y x &key (rcond -1))
;;   (lb-call ((common-type (common-float-type x y))
;;             (real-type (real-lla-type common-type))
;;             (procedure (lb-procedure-name common-type gelsd))
;;             ((:matrix x% common-type (m m%) (n n%) :output :copy) x)
;;             ((:matrix y% common-type m2 (nrhs nrhs%) :output b) y)
;;             ;; !!! fix for IWORK in DGELSD, need to remove once it is fixed in
;;             ;; !!! LAPACK.  Currently using a conservative estimate, as if
;;             ;; !!! SMLSIZ=0.
;;             (minmn (min m n))
;;             (nlvl (max 0 (1+ (ceiling (log minmn 2)))))
;;             ((:work iwork% :integer) (* minmn (+ 11 (* 3 nlvl))))
;;             ;; !!! end of fix, comment out iwork% below when removed.
;;             ((:work-queries lwork% 
;;                             (work% common-type)
;;                             (rwork% real-type t)
;;                             ;; (iwork% :integer)
;;                             ))
;;             ((:check info%))
;;             ((:output s% real-type s) (min m n))
;;             ((:output rank% :integer rank) 1)
;;             ((:atom rcond% real-type) (coerce* rcond real-type)))
;;     (assert (= m m2))
;;     (unless (<= n m)
;;       (error 'not-enough-columns))
;;     (if (lla-complex? common-type)
;;         (call procedure m% n% nrhs% x% m% y% m% s% rcond% rank% work% lwork% 
;;               rwork% iwork% info%)
;;         (call procedure m% n% nrhs% x% m% y% m% s% rcond% rank% work% lwork% 
;;               iwork% info%))
;;     (values 
;;       (matrix-from-first-rows b n nrhs m)
;;       (sum-last-rows b m nrhs n)
;;       (- m n)
;;       `(:s ,s :rank ,(aref rank 0)))))

;;; !! maybe a generic function, to calculate Cholesky decompositions from SVD-based
;;; !! least-squares too

(defgeneric invert-xx (xx)
  (:documentation "Calculate (X^T X)-1 (which is used for calculating the
variance of estimates) and return as a decomposition.  Usually XX is a
decomposition itself, eg QR returned by least squares.  Note: this can be used
to generate random draws, etc."))

(defmethod invert-xx ((qr qr))
  ;; Notes: X = QR, and thus X^T X = R^T Q^T Q R = R^T R because Q is
  ;; orthogonal, also (X^T X)^-1 = R^-1 (R^T)-1
  (bind (((:slots r) qr))
    (assert (<= (ncol r) (nrow r)))
    (make-instance 'matrix-square-root :left-square-root (invert r))))

;;;; constrained-least-squares

;; (defgeneric constrained-least-squares (y x z w)
;;   (:documentation "Solve the (linearly) constrained least squares
;;   problem min_b |y-Xb|_2 subject to Zx=w.  Return b."))

;; (defmethod constrained-least-squares ((y vector) (x dense-matrix-like)
;;                                       (z dense-matrix-like) (w vector))
;;   "Solve the (linearly) constrained least squares problem min_b |y-Xb|_2
;;   subject to Zx=w."
;;   ;; Note: mapping between the function parameters/variables and
;;   ;; LAPACK counterparts is as follows: y->c, X->A, x->b, Z->B, w->d,
;;   (lb-call ((common-type (common-float-type y x z w))
;;             (procedure (lb-procedure-name common-type gglse))
;;             ((:matrix x% common-type (m m%) (n n%) :output :copy) x)
;;             ((:matrix z% common-type (p p%) n2 :output :copy) z)
;;             ((:vector y% common-type m2 :copy) y)
;;             ((:vector w% common-type p2 :copy) w)
;;             ((:output b% common-type b) n)
;;             ((:work-queries lwork% (work% common-type)))
;;             ((:check info%)))
;;     ;; !! after call, x and z contain decompositions, currently not collected
;;     ;; !! after call, y contains sum of squares, currently not collected
;;     (assert (= n n2) () "Dimension mismatch between z and x")
;;     (assert (= m m2) () "Dimension mismatch between x and y")
;;     (assert (= p p2) () "Dimension mismatch between z and w")
;;     (call procedure m% n% p% x% m% z% p% y% w% b% work% lwork% info%)
;;     b))


;; ;;;; Cholesky factorization

(defgeneric cholesky (a)
  (:documentation "Cholesky factorization."))

(eval-when (:compile-toplevel :load-toplevel :execute)
  (assert (eql (load-time-value (hermitian-orientation :lapack)) +l+)
          () 
          "Cholesky implementation below assumes that hermitian orientation is
          L."))

(defmethod cholesky ((a hermitian-matrix))
  (lb-call ((a (elements a))
            (common-type (common-float-type a))
            ((:lapack potrf common-type))
            ((:array a% common-type :dimensions (n n2) :output l) a))
    (assert (= n n2))
    (call +l+ n a% n)
    (make-instance 'cholesky :left-square-root 
                   (make-matrix :lower nil :initial-contents l))))

(defmethod left-square-root ((hermitian-matrix hermitian-matrix))
  (left-square-root (cholesky hermitian-matrix)))

;; ;;; SVD

;; (defgeneric svd (a &key left right)
;;   (:documentation "Return singular value decomposition A, with left- and right
;; singular vectors as second and third values (when requested, see vector
;; specifications).  Valid vector specifications are :NONE, :SINGULAR (singular
;; vectors only) and :ALL.  Return values S (singular values, descending order, as
;; a DIAGONAL), U (left singular vectors, DENSE-MATRIX, or NIL), VT ([conjugate]
;; transpose of right singular vectors, DENSE-MATRIX, or NIL)."))

;; (defmethod svd ((a dense-matrix-like) &key (left :none) (right :none))
;;   (lb-call ((type (common-float-type a))
;;             (real-type (real-lla-type type))
;;             (complex? (lla-complex? type))
;;             (procedure (lb-procedure-name type gesvd))
;;             ((:matrix a% type (m m%) (n n%) :output :copy) a)
;;             (min-mn (min m n))
;;             ((:values u-ncol jobu) (ecase left
;;                                      (:none (values 1 #\N)) ; LAPACK needs >=1
;;                                      (:singular (values min-mn #\S))
;;                                      (:all (values m #\A))))
;;             ((:values vt-nrow jobvt) (ecase right
;;                                        (:none (values 1 #\N)) ; LAPACK needs >=1
;;                                        (:singular (values min-mn #\S))
;;                                        (:all (values n #\A))))
;;             ((:char jobu%) jobu)
;;             ((:char jobvt%) jobvt)
;;             ((:integer vt-nrow%) vt-nrow)
;;             ((:output s% real-type s) min-mn)
;;             ((:output u% type u) (* m u-ncol))
;;             ((:output vt% type vt) (* vt-nrow n))
;;             ((:work rwork% real-type) (if complex? (* 5 min-mn) 0))
;;             ((:work-queries lwork% (work% type)))
;;             ((:check info%)))
;;     (with-lapack-traps-masked
;;       (if complex?
;;           (call procedure jobu% jobvt% m% n% a% m% s% u% m% vt% vt-nrow%
;;                 work% lwork% rwork% info%)
;;           (call procedure jobu% jobvt% m% n% a% m% s% u% m% vt% vt-nrow%
;;                 work% lwork% info%)))
;;     (values (make-diagonal% s)
;;             (unless (eq left :none) (make-matrix% m u-ncol u))
;;             (unless (eq right :none) (make-matrix% vt-nrow n vt)))))

;;; trace

;; (defgeneric tr (a)
;;   (:documentation "Trace of a square matrix.")
;;   (:method ((a dense-matrix-like))
;;     (check-type a square-matrix)
;;     (bind (((:lla-matrix a) a))
;;       (iter
;;         (for i :from 0 :below (nrow a))
;;         (summing (a (a-index i i))))))
;;   (:method ((a diagonal))
;;     (reduce #'+ (elements a))))

;; rank

;; (defun rank (matrix &key threshold (logrc-threshold 0.5))
;;   "Calculate the rank of a matrix using singular values - only those
;; above THRESHOLD in absolute value are counted.  If NIL, THRESHOLD is
;; determined automatically as (* (MAX NROW NCOL) FLOAT-EPSILON) times
;; the largest absolute singular value.  If (ABS (LOG (/ NROW NCOL))) is
;; above LOGRC-THRESHOLD, use (MM MATRIX T) or its transpose for faster
;; calculations; this can be disabled by setting LOG-THRESHOLD to NIL.
;; Returns the absolute values of the singular values as the second
;; value."
;;   (check-type logrc-threshold (or null (and number (satisfies plusp))))
;;   (bind (((:accessors-r/o nrow ncol) matrix)
;;          (ratio (log (/ nrow ncol)))
;;          ((:values matrix squared?)
;;           (cond 
;;             ((aand logrc-threshold (< ratio (- it))) ; fewer rows than columns
;;              (values (mm matrix t) t))
;;             ((aand logrc-threshold (< it ratio)) ; more rows than columns
;;              (values (mm t matrix) t))
;;             (t                          ; no matrix multiplication
;;              (values matrix nil))))
;;          (d (svd matrix))
;;          (d-elements (elements d)))
;;     (map-into d-elements (if squared? (lambda (x) (sqrt (abs x))) #'abs)
;;               d-elements)
;;     (let ((threshold (aif threshold
;;                           it
;;                           (* (max nrow ncol) (epsilon* (common-float-type d))
;;                              (reduce #'max d-elements)))))
;;       (values (count-if (lambda (x) (<= threshold (abs x))) (elements d)) d))))

;;; determinants

(defgeneric logdet (matrix)
  (:documentation "Logarithm of the determinant of a matrix.  Return -1, 1 or
  0 (or equivalent) to correct for the sign, as a second value."))

(defun det (matrix)
  "Determinant of a matrix.  If you need the log of this, use LOGDET
  directly."
  (bind (((:values logdet sign) (logdet matrix)))
    (if (zerop sign)
        0
        (* sign (exp logdet)))))

(defmacro log-with-sign% (value sign-changes block-name)
  "Log of (ABS VALUE), increments SIGN-CHANGES when negative, return-from
block-name (values nil 0) when zero."
  (once-only (value)
    `(log (cond
            ((zerop ,value) (return-from ,block-name (values nil 0)))
            ((minusp ,value) (incf ,sign-changes) (- ,value))
            (t ,value)))))

(defun diagonal-log-sum% (matrix &optional (sign-changes 0))
  "Sum of the log of the elements in the diagonal.  Sign-changes counts the
negative values, and may be started at something else than 0 (eg in case of
pivoting).  Return (values NIL 0) in case of encountering a 0."
  (bind (((nrow ncol) (array-dimensions matrix)))
    (assert (= nrow ncol))
    (values 
      (iter
        (for i :from 0 :below nrow)
        (summing (log-with-sign% (aref matrix i i)
                                 sign-changes diagonal-log-sum%)))
      (if (evenp sign-changes) 1 -1))))

(defmethod logdet ((matrix array))
  (let* ((lu (lu matrix)))
    (diagonal-log-sum% (lu lu) (permutations lu))))

(defmethod logdet ((matrix lower-triangular-matrix))
  (diagonal-log-sum% (elements matrix)))

(defmethod logdet ((matrix upper-triangular-matrix))
  (diagonal-log-sum% (elements matrix)))

;; (defmethod logdet ((matrix hermitian-matrix))
;;   (let ((sign-changes 0))
;;     (values 
;;       (reduce #'+ (eigen matrix)
;;               :key (lambda (e)
;;                      (log-with-sign% e sign-changes logdet)))
;;       (if (evenp sign-changes) 1 -1))))

;; (defun matrix-cond (matrix)
;;   "Calculate the condition number of a matrix (with respect to the Euclidean
;;   norm).  Implemented as the ratio of the largest and smallest singular values."
;;   (check-type matrix dense-matrix-like)
;;   (let* ((s (elements (svd matrix)))
;;          (s-min (aref s (1- (length s)))))
;;     (if (zerop s-min)
;;         (error "Matrix is not full rank.")
;;         (/ (aref s 0) s-min))))
