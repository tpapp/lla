(in-package :lla)

;;;; Higher level linear algebra functions defined here.
;;;
;;; General convention for vectors in places of matrices: should be interpreted as a conforming vector.  If the result is an 1xn or nx1 matrix, it should be converted to a vector iff some related argument was a vector.  For example, (solve a b) should be a vector iff b is a vector, otherwise it should be a matrix.
;;;
;;;;  matrix multiplication
;;;
;;; The general matrix multiplication function is MMM.  It can be called with more than two arguments: it multiplies matrices from the left, using MM.
;;;
;;; MM takes two matrices and an optional scalar.  Specialized methods may allow matrices to be replaced by symbols etc, in which case they have a special interpretation (eg T stands for "multiply matrix by conjugate transpose").
;;;
;;; !! Various optimizations, deferred to the future:
;;; - compiler macros could notice (transpose matrix)
;;; - mmm could detect and collect scalars
;;; - methods for triangular matrices

(defun dimensions-as-matrix (array orientation)
  "If ARRAY is a vector, return its dimensions as a matrix with given ORIENTATION (:ROW or :COLUMN) as multiple values, and T as the third value.  If it is matrix, just return the dimensions and NIL.  Used for considering vectors as conforming matrices (eg see MM)."
  (ecase (array-rank array)
    (1 (let ((d0 (array-dimension array 0)))
         (ecase orientation
           (:row (values 1 d0 t))
           (:column (values d0 1 t)))))
    (2 (let+ (((d0 d1) (array-dimensions array)))
         (values d0 d1 nil)))))

(defgeneric mm (a b)
  (:documentation  "Matrix multiplication of A and B.")
  (:method (a b) (mm (aops:as-array a) (aops:as-array b))))

(defmethod aops:as-array ((matrix-square-root matrix-square-root))
  (mm (left-square-root matrix-square-root) t))

(defun mmm (&rest matrices)
  "Multiply arguments from left to right using MM."
  (reduce #'mm matrices))

(defmethod mm ((a array) (b array))
  (let+ ((common-type (common-float-type a b))
         ((&values a0 a1 av) (dimensions-as-matrix a :row))
         ((&values b0 b1 bv) (dimensions-as-matrix b :column))
         (c-dimensions (cond
                         ((and av bv) (error "MM called with two vectors."))
                         ((or av bv) (* a0 b1))
                         (t (list a0 b1)))))
    (assert (= a1 b0) () "Incompatible dimensions.")
    ;; here C=AB <=> C^T=B^T A^T, so in the argument list, A and B are
    ;; interchanged
    (blas-call ("gemm" common-type c)
      #\N #\N (&integers b1 a0 b0) 1 (&array-in b) (&integer b1)
      (&array-in a) (&integer a1) 0
      (&array-out (&new c) :dimensions c-dimensions :type common-type)
      (&integer b1))))

;;; !! this is how we could speed things up with compiler macros: have a
;;; !! compiler macro transform (MM (TRANSPOSE FOO) BAR) to (MM-TN FOO BAR),
;;; !! which is a generic function and would use MM-INTERNAL accordingly.  Or,
;;; !! we could have (lazy-transpose* FOO) explicitly, and have MM dispatch on
;;; !! that.  Will investigate.

(defun mm-hermitian% (a transpose-left?)
  "Calculate A*op(A) if TRANSPOSE-LEFT? is NIL, and op(A)*A otherwise.  op() is always conjugate transpose, but may be implemented as a transpose if A is real, in which case the two are equivalent.  This function is meant to be used internally, and is not exported."
  ;; implementation note: no transpose is necessary
  (let+ (((a0 a1) (array-dimensions a))
         ((&values dim-c other-dim-a)
          (if transpose-left?
              (values a1 a0)
              (values a0 a1)))
         (type (common-float-type a)))
    (blas-call (("syrk" "herk") type (hermitian-matrix c))
      #\U (&char (if transpose-left? #\N #\C))
      (&integers dim-c other-dim-a) 1
      (&array-in a) (&integer a1) 0
      (&array-out (&new c) :dimensions (list dim-c dim-c) :type type)
      (&integer dim-c))))

(defmethod mm ((a array) (b (eql t)))
  ;; A A^T
  (mm-hermitian% a nil))

(defmethod mm ((a (eql t)) (b array))
  ;; A^T A
  (mm-hermitian% b t))

(defmethod mm ((a wrapped-matrix) (b wrapped-matrix))
  (mm (aops:as-array a) (aops:as-array b)))

(defmethod mm ((a wrapped-matrix) b)
  (mm (aops:as-array a) b))

(defmethod mm (a (b wrapped-matrix))
  (mm a (aops:as-array b)))

;;; (mm vector t) is the dot product

(defmethod mm ((a vector) (b vector))
  (assert (= (length a) (length b)))
  (loop for a-elt :across a
        for b-elt :across b
        summing (* a-elt (conjugate b-elt))))

(defmethod mm ((a vector) (b (eql t)))
  (declare (inline absolute-square))
  (reduce #'+ a :key #'absolute-square))

(defmethod mm ((a (eql t)) (b vector))
  (declare (inline absolute-square))
  (reduce #'+ b :key #'absolute-square))

;;; outer product

(defgeneric outer (a b)
  (:documentation "Return the outer product column(a) row(b)^H.  If either A or B is T, they are taken to be conjugate transposes of the other argument.")
  (:method ((a vector) (b (eql t)))
    (mm (aops:reshape-col a) t))
  (:method ((a (eql t)) (b vector))
    (mm (aops:reshape-col b) t))
  (:method ((a vector) (b vector))
    (let* ((a0 (length a))
           (b0 (length b))
           ;; !!! currently we are not using the narrowest element type
           (result (make-array (list a0 b0)))
           (index 0))
      (dotimes (a-index a0)
        (let ((a-elt (aref a a-index)))
          (dotimes (b-index b0)
            (setf (row-major-aref result index)
                  (* a-elt (conjugate (aref b b-index))))
            (incf index))))
      result)))

;;; mm for diagonal matrices

(defmethod mm ((a diagonal-matrix) (b array))
  (let+ (((b0 b1) (array-dimensions b))
         (a (diagonal-matrix-elements a))
         ;; !!! currently we are not using the narrowest element type
         (c (make-array (list b0 b1)
                        :element-type (elementwise-float-contagion a b)))
         (i 0))
    (assert (= (length a) b0) () 'lla-incompatible-dimensions)
    (dotimes (row b0)
      (let ((a (row-major-aref a row)))
        (dotimes (col b1)
          (setf (row-major-aref c i) (* (row-major-aref b i) a))
          (incf i))))
    c))

(defmethod mm ((a array) (b diagonal-matrix))
  (let+ (((a0 a1) (array-dimensions a))
         (b (diagonal-matrix-elements b))
         ;; !!! currently we are not using the narrowest element type
         (c (make-array (list a0 a1)
                        :element-type (elementwise-float-contagion a b)))
         (i 0))
    (assert (= a1 (length b)) () 'lla-incompatible-dimensions)
    (dotimes (row a0)
      (dotimes (col a1)
        (setf (row-major-aref c i)
              (* (row-major-aref a i) (row-major-aref b col)))
        (incf i)))
    c))

(defmethod mm ((a vector) (b diagonal-matrix))
  (e* a (diagonal-matrix-elements b)))

(defmethod mm ((a diagonal-matrix) (b vector))
  (e* (diagonal-matrix-elements a) b))

(defmethod mm ((a diagonal-matrix) (b diagonal-matrix))
  (diagonal-matrix (e* (diagonal-matrix-elements a) (diagonal-matrix-elements b))))

(defmethod mm ((a diagonal-matrix) (b (eql t)))
  ;; !!! currently we are not using the narrowest element type
  (let+ (((&accessors-r/o diagonal-matrix-elements) a))
    (diagonal-matrix (map (type-of diagonal-matrix-elements)
                          #'absolute-square diagonal-matrix-elements))))

(defmethod mm ((a (eql t)) (b diagonal-matrix))
  (mm b a))


;;; hermitian (symmetric) updates

;; (defun update-hermitian (a x &optional (alpha 1))
;;   "Return alpha x x^H + A.  A has to be a hermitian matrix, and ALPHA real."
;;   (check-type a hermitian-matrix)
;;   (check-type x vector)
;;   (lb-call ((common-type (lb-target-type a x alpha))
;;             (real-type (absolute-square-type common-type))
;;             (procedure (lb-procedure-name common-type syr her))
;;             ((:matrix a% common-type (n n%) n2 :output c) a)
;;             ((:vector x% common-type n3) x)
;;             ((:atom alpha% real-type) (coerce* alpha real-type))
;;             ((:char u%) #\U)
;;             ((+integer+ one%) 1))
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
;;             ((+integer+ one%) 1))
;;     (assert (= n n2 n3 n4))
;;     (call procedure u% n% alpha% x% one% y% one% a% n%)
;;     (make-matrix% n n c :kind :hermitian)))


;;;; LU factorization

(defgeneric lu (a)
  (:documentation "LU decomposition of A"))

(defmethod lu ((a array))
  (let+ (((a0 a1) (array-dimensions a)))
    (lapack-call ("getrf" (common-float-type a)
                          (make-instance 'lu :lu lu :ipiv ipiv))
      (&integers a0 a1)
      (&array-in/out
          (:input a :transpose? t)
          (:output (&new lu) :transpose? t))
      (&integer a0) (&array-out (&new ipiv) :dimensions (min a0 a1) :type +integer+)
      (&info :condition nil))))

;;;; Hermitian factorization

(defgeneric hermitian-factorization (a)
  (:documentation "Compute the hermitian factorization."))

(defmethod hermitian-factorization ((a hermitian-matrix))
  (let+ ((a (wrapped-matrix-elements a))
         ((a0 a1) (array-dimensions a)))
    (assert (= a0 a1) () 'lla-incompatible-dimensions)
    (lapack-call-w/query (("sytrf" "hetrf") (common-float-type a)
                          (make-instance 'hermitian-factorization
                                         :factor factor :ipiv ipiv))
      #\U
      (&integers a0)
      (&array-in/out (:input a) (:output (&new factor)))
      (&integer a0)
      (&array-out (&new ipiv) :dimensions a0 :type +integer+)
      (&work-query) (&info))))

;;;; solving linear equations

(defgeneric solve (a b)
  (:documentation "Return X that solves AX=B.  When B is a vector, so is X.")
  (:argument-precedence-order b a))

(defmethod solve (a (b wrapped-matrix))
  (solve a (aops:as-array b)))

(defmethod solve ((lu lu) (b array))
  (let+ (((&slots lu ipiv) lu)
         ((lu0 lu1) (array-dimensions lu))
         ((&values b0 b1 &ign) (dimensions-as-matrix b :column)))
    (assert (= lu0 lu1 b0) () 'lla-incompatible-dimensions)
    (lapack-call ("getrs" (common-float-type lu b) x)
      #\N (&integer lu0) (&integer b1)
      (&array-in lu :transpose? t)
      (&integer lu0) (&array-in ipiv :type +integer+)
      (&array-in/out
          (:input b :transpose? t)
          (:output (&new x) :transpose? t))
      (&integer lu0) &info)))

(defmethod solve ((a array) (b array))
  (let+ (((a0 a1) (array-dimensions a))
         ((&values b0 b1 &ign) (dimensions-as-matrix b :column)))
    (assert (= a0 a1 b0) () 'lla-incompatible-dimensions)
    (lapack-call ("gesv" (common-float-type a b) x)
      (&integer a0) (&integer b1) (&array-in a :transpose? t)
      (&integer a0)
      (&work a0 +integer+)
      (&array-in/out
          (:input b :transpose? t)
          (:output (&new x) :transpose? t))
      (&integer a0) &info)))

(defmethod solve ((cholesky cholesky) b)
  (let+ ((a (wrapped-matrix-elements (left-square-root cholesky)))
         ((&values b0 b1 &ign) (dimensions-as-matrix b :column))
         ((a0 a1) (array-dimensions a)))
    (assert (= a0 a1 b0) () 'lla-incompatible-dimensions)
    (lapack-call ("potrs" (common-float-type a b) x)
      #\U (&integers a0 b1) (&array-in a) (&integer a0)
      (&array-in/out
          (:input b :transpose? t)
          (:output (&new x) :transpose? t))
      (&integer b0) &info)))

(defmethod solve ((hermitian-matrix hermitian-matrix) b)
  (let+ ((a (wrapped-matrix-elements hermitian-matrix))
         ((&values b0 b1) (dimensions-as-matrix b :column))
         ((a0 a1) (array-dimensions a)))
    (assert (= a0 a1 b0) () 'lla-incompatible-dimensions)
    (lapack-call ("posv" (common-float-type a b) x)
      #\U (&integers a0 b1)
      (&array-in a :force-copy? t)
      (&integer a0)
      (&array-in/out
          (:input b :transpose? t)
          (:output (&new x) :transpose? t))
      (&integer b0) &info)))

(defmethod solve ((a hermitian-factorization) (b array))
  (let+ (((&slots-r/o factor ipiv) a)
         ((a0 a1) (array-dimensions factor))
         ((&values b0 b1) (dimensions-as-matrix b :column)))
    (assert (= a0 a1 b0) () 'lla-incompatible-dimensions)
    (lapack-call (("sytrs" "hetrs") (common-float-type factor b) x)
      #\U (&integers a0 b1) (&array-in factor) (&integer a0)
      (&array-in ipiv :type +integer+)
      (&array-in/out
          (:input b :transpose? t)
          (:output (&new x) :transpose? t))
      (&integer b0) &info)))

(defmethod solve ((a hermitian-matrix) (b array))
  (solve (hermitian-factorization a) b))

(defun trsm% (a a-upper? b)
  "Wrapper for BLAS routine xTRSM.  Solve AX=B, where A is triangular."
  (let+ (((a0 a1) (array-dimensions a))
         ((&values b0 b1 &ign) (dimensions-as-matrix b :column)))
    (assert (= a0 a1 b0) () 'lla-incompatible-dimensions)
    (blas-call ("trsm" (common-float-type a b) x)
      #\R (&char (if a-upper? #\L #\U)) #\N #\N (&integers b1 b0) 1
      (&array-in a) (&integer a1) (&array-in/out (:input b) (:output (&new x)))
      (&integer b1))))

(defmethod solve ((a lower-triangular-matrix) b)
  (trsm% (wrapped-matrix-elements a) nil b))

(defmethod solve ((a upper-triangular-matrix) b)
  (trsm% (wrapped-matrix-elements a) t b))

(defmethod solve ((a diagonal-matrix) b)
  (mm (e/ a) b))

(defgeneric invert (a &key &allow-other-keys)
  (:documentation "Invert A.  The inverse of matrix factorizations are other factorizations when appropriate, otherwise the result is a matrix.  Usage note: inverting dense matrices is unnecessary and unwise in most cases, because it is numerically unstable.  If you are solving many Ax=b equations with the same A, use a matrix factorization like LU."))

(defmethod invert ((a array) &key)
  (invert (lu a)))

(defmethod invert ((lu lu) &key)
  (let+ (((&slots-r/o lu ipiv) lu)
         ((lu0 lu1) (array-dimensions lu)))
    (assert (= lu0 lu1 (length ipiv)))
    (lapack-call-w/query ("getri" (common-float-type lu) inverse)
      (&integer lu0)
      (&array-in/out (:input lu :transpose? t)
          (:output (&new inverse) :transpose? t))
      (&integer lu0) (&array-in ipiv :type +integer+) (&work-query) &info)))

(defmethod invert ((a hermitian-factorization) &key)
  (let+ (((&slots-r/o factor ipiv) a)
         ((a0 a1) (array-dimensions factor)))
    (assert (= a0 a1))
    (lapack-call (("sytri" "hetri") (common-float-type factor)
                  (hermitian-matrix inverse))
      #\U (&integer a0)
      (&array-in/out (:input factor) (:output (&new inverse)))
      (&integer a0)
      (&array-in ipiv :type +integer+) (&work a0) &info)))

(defmethod invert ((a hermitian-matrix) &key)
  (invert (hermitian-factorization a)))

(defun invert-triangular% (a upper? unit-diag?)
  "Invert a dense (triangular) matrix using the LAPACK routine *TRTRI.  UPPER? indicates if the matrix is in the upper or the lower triangle of a matrix, UNIT-DIAG? indicates whether the diagonal is supposed to consist of 1s.  For internal use, not exported."
  (let+ (((a0 a1) (array-dimensions a)))
    (assert (= a0 a1))
    (lapack-call ("trtri" (common-float-type a) inverse)
      (&char (if upper? #\L #\U))
      (&char (if unit-diag? #\U #\N))
      (&integer a0)
      (&array-in/out (:input a) (:output (&new inverse)))
      (&integer a0)
      &info)))

(defmethod invert ((a upper-triangular-matrix) &key)
  (upper-triangular-matrix
   (invert-triangular% (wrapped-matrix-elements a) t nil)))

(defmethod invert ((a lower-triangular-matrix) &key)
  (lower-triangular-matrix
   (invert-triangular% (wrapped-matrix-elements a) nil nil)))

(defmethod invert ((a cholesky) &key)
  (let+ ((a (wrapped-matrix-elements (left-square-root a)))
         ((a0 a1) (array-dimensions a)))
    (assert (= a0 a1))
    (lapack-call ("potri" (common-float-type a)
                          (hermitian-matrix inverse))
      #\U (&integer a0) (&array-in/out (:input a) (:output (&new inverse)))
      (&integer a0) &info)))

(defmethod invert ((diagonal diagonal-matrix) &key (tolerance 0))
  "For pseudoinverse, suppressing diagonal elements below TOLERANCE \(if given, otherwise / is just used without any checking."
  (let ((elements (diagonal-matrix-elements diagonal)))
    (diagonal-matrix (map `(simple-array
                            ,(clnu::elementwise-float-contagion elements) (*))
                          (cond
                            ((null tolerance) #'/)
                            ((and (numberp tolerance) (<= 0 tolerance))
                             (lambda (x)
                               (if (<= (abs x) tolerance) 0 (/ x))))
                            ((and (numberp tolerance) (zerop tolerance))
                             (lambda (x)
                               (if (zerop x) 0 (/ x))))
                            (t (error "Invalid tolerance argument.")))
                          elements))))

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
;;             (real-type (absolute-square-type common-type))
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

(defun last-rows-ss (matrix nrhs common-type)
  "Calculate the sum of squares of the last rows of MATRIX columnwise,
omitting the first NRHS rows.  If MATRIX is a vector, just do this for the last elements.  Used for interfacing with xGELS."
  (ecase (array-rank matrix)
    (1 (reduce #'+ matrix :key #'absolute-square :start nrhs))
    (2 (let+ (((m n) (array-dimensions matrix))
              (element-type (lisp-type (absolute-square-type common-type)))
              (sum (make-array n :element-type element-type
                               :initial-element (coerce 0 element-type)))
              (matrix-index (array-row-major-index matrix nrhs 0)))
         (loop repeat (- m nrhs) do
           (dotimes (sum-index n)
             (incf (aref sum sum-index)
                   (absolute-square (row-major-aref matrix matrix-index)))
             (incf matrix-index)))
         sum))))

(defun least-squares-qr (y x &key &allow-other-keys)
  "Least squares using QR decomposition.  Additional values returned: the QR decomposition of X.  See LEAST-SQUARES for additional documentation.  Usage note: SVD-based methods are recommended over this one, unless X is well-conditioned."
  ;; Note: the naming convention (y,X,b) is different from LAPACK's
  ;; (b,A,x).  Sorry if this creates confusion, I decided to follow
  ;; standard statistical notation.
  ;;
  ;; implementation note: there is no point in transposing the last rows of
  ;; b-and-ss, we could just sum them as is, would even be better for
  ;; linearity of memory access
  (let+ (((&values y0 y1) (dimensions-as-matrix y :column))
         ((x0 x1) (array-dimensions x))
         (df (- x0 x1))
         (common-type (common-float-type y x)))
    (assert (= y0 x0) () 'lla-incompatible-dimensions)
    (assert (plusp df) () 'not-enough-columns)
    (lapack-call-w/query
        ("gels"
         common-type
         (values
          (copy-array (aops:partition b-and-ss 0 x1))
          (last-rows-ss b-and-ss x1 common-type)
          df
          (make-instance 'qr :qr qr)))
      #\N (&integers x0 x1 y1)
      (&array-in/out
          (:input x :transpose? t)
          (:output (&new qr) :transpose? t))
      (&integer x0)
      (&array-in/out
          (:input y :transpose? t)
          (:output (&new b-and-ss) :transpose? t))
      (&integer y0) (&work-query) &info)))

(defmethod qr ((a array))
  (let+ (((a0 a1) (array-dimensions a)))
    (lapack-call-w/query
        ("geqrf" (common-float-type a)
                 (make-instance 'qr :qr qr :tau tau))
      (&integers a0 a1)
      (&array-in/out
          (:input a :transpose? t)
          (:output (&new qr) :transpose? t))
      (&integer a0) (&array-out (&new tau) :dimensions (min a0 a1)) (&work-query)
      &info)))

;; (defun least-squares-svd-d (y x &key (rcond -1))
;;   (lb-call ((common-type (common-float-type x y))
;;             (real-type (absolute-square-type common-type))
;;             (procedure (lb-procedure-name common-type gelsd))
;;             ((:matrix x% common-type (m m%) (n n%) :output :copy) x)
;;             ((:matrix y% common-type m2 (nrhs nrhs%) :output b) y)
;;             ;; !!! fix for IWORK in DGELSD, need to remove once it is fixed in
;;             ;; !!! LAPACK.  Currently using a conservative estimate, as if
;;             ;; !!! SMLSIZ=0.
;;             (minmn (min m n))
;;             (nlvl (max 0 (1+ (ceiling (log minmn 2)))))
;;             ((:work iwork% +integer+) (* minmn (+ 11 (* 3 nlvl))))
;;             ;; !!! end of fix, comment out iwork% below when removed.
;;             ((:work-queries lwork%
;;                             (work% common-type)
;;                             (rwork% real-type t)
;;                             ;; (iwork% +integer+)
;;                             ))
;;             ((:check info%))
;;             ((:output s% real-type s) (min m n))
;;             ((:output rank% +integer+ rank) 1)
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
  (:documentation "Calculate (X^T X)-1 (which is used for calculating the variance of estimates) and return as a decomposition.  Usually XX is a decomposition itself, eg QR returned by least squares.  Note: this can be used to generate random draws, etc."))

(defmethod invert-xx ((qr qr))
  ;; Notes: X = QR, and thus X^T X = R^T Q^T Q R = R^T R because Q is
  ;; orthogonal, also (X^T X)^-1 = R^-1 (R^T)-1
  (let+ ((r (qr-r qr)))
    (assert (<= (aops:ncol r) (aops:nrow r)))
    (xx (invert r))))

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


;;;; Cholesky factorization

(defgeneric cholesky (a)
  (:documentation "Cholesky factorization."))

(defmethod cholesky ((a hermitian-matrix))
  (let+ ((a (wrapped-matrix-elements a))
         ((a0 a1) (array-dimensions a)))
    (assert (= a0 a1))
    (lapack-call ("potrf" (common-float-type a)
                          (make-cholesky (lower-triangular-matrix l)))
      #\U (&integer a0)
      (&array-in/out (:input a) (:output (&new l)))
      (&integer a0) &info)))

(defmethod left-square-root ((hermitian-matrix hermitian-matrix))
  (left-square-root (cholesky hermitian-matrix)))

(defmethod left-square-root ((a diagonal-matrix))
  (esqrt a))

;;; spectral factorization

(defun eigenvalues (a &key (abstol 0))
  "Return the eigenvalues of A.  See the documentation of
SPECTRAL-FACTORIZATION about ABSTOL."
  (check-type a hermitian-matrix)
  (let+ ((a (wrapped-matrix-elements a))
         (type (common-float-type a))
         (real-type (absolute-square-type type))
         ((a0 a1) (array-dimensions a)))
    (assert (= a0 a1))
    (with-fp-traps-masked
      (if (complex? type)
          (error "needs to be written, report this as an issue")
          (lapack-call-w/query (("syevr" "heevr") type w)
            #\N #\A #\U (&integer a0)
            (&array-in/out (:input a) ())
            (&integer a0) nil nil nil nil (&atom abstol :type real-type)
            (&work 1 +integer+) (&array-out (&new w) :dimensions a0 :type real-type)
            nil (&integer a0) (&work (* 2 (max 1 a0))) (&work-query)
            (&work-query +integer+) &info)))))

(defun spectral-factorization (a &key (abstol 0))
  "Return a spectral factorization of A.

The LAPACK manual says the following about ABSTOL:

The absolute error tolerance for the eigenvalues.  An approximate eigenvalue is accepted as converged when it is determined to lie in an interval [a,b] of width less than or equal to

                  ABSTOL + EPS *   max( |a|,|b| ) ,

where EPS is the machine precision.  If ABSTOL is less than or equal to zero, then EPS*|T| will be used in its place, where |T| is the 1-norm of the tridiagonal matrix obtained by reducing A to tridiagonal form.

See \"Computing Small Singular Values of Bidiagonal Matrices with Guaranteed High Relative Accuracy,\" by Demmel and Kahan, LAPACK Working Note #3.

If high relative accuracy is important, set ABSTOL to DLAMCH( 'Safe minimum').  Doing so will guarantee that eigenvalues are computed to high relative accuracy when possible in future releases.  The current code does not make any guarantees about high relative accuracy, but furutre releases will. See J. Barlow and J. Demmel, \"Computing Accurate Eigensystems of Scaled Diagonally Dominant Matrices\", LAPACK Working Note #7, for a discussion of which matrices define their eigenvalues to high relative accuracy."
  (check-type a hermitian-matrix)
  (let+ ((a (wrapped-matrix-elements a))
         (type (common-float-type a))
         (real-type (absolute-square-type type))
         ((a0 a1) (array-dimensions a)))
    (assert (= a0 a1))
    (with-fp-traps-masked
      (if (complex? type)
          (error "needs to be written, report this as an issue")
          (lapack-call-w/query (("syevr" "heevr") type
                                (make-spectral-factorization
                                 :z z :w (diagonal-matrix w)))
            #\V #\A #\U (&integer a0) (&array-in/out (:input a) ())
            (&integer a0) nil nil nil nil (&atom abstol :type real-type)
            (&work 1 +integer+) (&array-out (&new w) :dimensions a0 :type real-type)
            (&array-out (&new z) :dimensions (list a0 a1)) (&integer a0)
            (&work (* 2 (max 1 a0))) (&work-query) (&work-query +integer+)
            &info)))))

(defmethod aops:as-array ((sf spectral-factorization))
  (let+ (((&structure-r/o spectral-factorization- z w) sf))
    (assert z)
    (mm (mm z (esqrt w)) t)))

;;; SVD

(defgeneric svd (a &optional vectors)
  (:documentation "Return singular value decomposition A.

  VECTORS determines how singular vectors are calculated:

  - NIL sets U and VT to NIL
  - :ALL makes U and VT square, with dimensions conforming to A
  - :THIN makes the larger of U and VT rectangular.  This means that not all
    of the singular vectors are calculates, and saves computational time when
    A is far from square."))

(defmethod svd ((a array) &optional (vectors :thin))
  (let+ (((a0 a1) (array-dimensions a))
         (min (min a0 a1))
         (type (common-float-type a)))
    (if (complex? type)
        ;; for the complex case, we have to transpose
        (error "not implemented yet, report as a bug")
        ;; for the real case, we don't have to transpose, just exchange order
        (let+ (((&values jobz u1 vt0) (ecase vectors
                                        ((nil) (values #\N 0 0))
                                        (:thin (values #\S min min))
                                        (:all (values #\A a0 a1)))))
          (lapack-call-w/query ("gesdd" type
                                        (make-svd :d (diagonal-matrix d)
                                                  :u (when vectors u)
                                                  :vt (when vectors vt)))
            (&char jobz)
            (&integers a1 a0)
            (&array-in a :force-copy? t)
            (&integer a1)
            (&array-out (&new d) :dimensions min)
            (&array-out (&new vt) :dimensions (list vt0 a1))
            (&integer (max a1 1))
            (&array-out (&new u) :dimensions (list a0 u1))
            (&integer (max u1 1))
            (&work-query) (&work (* 8 min) +integer+)
            &info)))))

(defmethod aops:as-array ((svd svd))
  (let+ (((&structure-r/o svd- u d vt) svd)
         (n (aops:nrow d)))
    (mmm (if (= (aops:ncol u) n)
             u
             (slice u t (cons 0 n)))
         d
         (if (= (aops:nrow vt) n)
             vt
             (slice vt (cons 0 n) t)))))

;;; trace

(defun sum-diagonal% (array)
  "Sum diagonal of array, checking that it is square."
  (let+ (((a0 a1) (array-dimensions array)))
    (assert (= a0 a1))
    (loop for i below a0 summing (aref array i i))))

(defgeneric tr (a)
  (:documentation "Trace of a square matrix.")
  (:method ((a array))
    (sum-diagonal% a))
  (:method ((a wrapped-matrix))
    (sum-diagonal% (wrapped-matrix-elements a)))
  (:method ((a diagonal-matrix))
    (sum (diagonal-matrix-elements a))))

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
;;          ((&values matrix squared?)
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
  (:documentation "Logarithm of the determinant of a matrix.  Return -1, 1 or 0 (or equivalent) to correct for the sign, as a second value."))

(defun det (matrix)
  "Determinant of a matrix.  If you need the log of this, use LOGDET
  directly."
  (let+ (((&values logdet sign) (logdet matrix)))
    (if (zerop sign)
        0
        (* sign (exp logdet)))))

(defmacro log-with-sign% (value sign-changes block-name)
  "Log of (ABS VALUE), increments SIGN-CHANGES when negative, return-from block-name (values nil 0) when zero."
  (once-only (value)
    `(log (cond
            ((zerop ,value) (return-from ,block-name (values nil 0)))
            ((minusp ,value) (incf ,sign-changes) (- ,value))
            (t ,value)))))

(defun diagonal-log-sum% (matrix &optional (sign-changes 0))
  "Sum of the log of the elements in the diagonal.  Sign-changes counts the negative values, and may be started at something else than 0 (eg in case of pivoting).  Return (values NIL 0) in case of encountering a 0."
  (let+ (((nrow ncol) (array-dimensions matrix)))
    (assert (= nrow ncol))
    (values
      (loop
        for i from 0 below nrow
        summing (log-with-sign% (aref matrix i i)
                                sign-changes diagonal-log-sum%))
      (if (evenp sign-changes) 1 -1))))

(defmethod logdet ((matrix array))
  (let* ((lu (lu matrix)))
    (diagonal-log-sum% (lu-matrix lu) (permutations lu))))

(defmethod logdet ((matrix lower-triangular-matrix))
  (diagonal-log-sum% (wrapped-matrix-elements matrix)))

(defmethod logdet ((matrix upper-triangular-matrix))
  (diagonal-log-sum% (wrapped-matrix-elements matrix)))

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
;;   (let* ((s (wrapped-matrix-elements (svd matrix)))
;;          (s-min (aref s (1- (length s)))))
;;     (if (zerop s-min)
;;         (error "Matrix is not full rank.")
;;         (/ (aref s 0) s-min))))
