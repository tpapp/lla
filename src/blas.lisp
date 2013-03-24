(in-package :lla)

;;;; Wrappers for BLAS linear algebra functions defined here.
;;;
;;; Only a select few BLAS functions are implemented currently. All
;;; functions take a SAFE? argument that defaults to *SAFE?* (T by
;;; default). With SAFE? true, arguments are checked to catch errors
;;; that BLAS would also catch and kill the lisp.

(defvar *safe?* t)

(defun gemm! (alpha a b beta c &key transpose-a? transpose-b?
              m n k lda ldb ldc (safe? *safe?*))
  "Basically C = ALPHA * A' * B' + BETA * C. A' is A or its transpose depending on TRANSPOSE-A?. B' is B or its transpose depending on TRANSPOSE-B?. Returns C.

A' is an MxK matrix. B' is a KxN matrix. C is an MxN matrix.

LDA is the width of the matrix A (not of A'). If A is not transposed, then K <= LDA, if it's transposed then M <= LDA.

LDB is the width of the matrix B (not of B'). If B is not transposed, then N <= LDB, if it's transposed then K <= LDB.

In the example below M=3, N=2, K=5, LDA=6, LDB=3, LDC=4. The cells marked with + do not feature in the calculation.

           N
          --+
          --+
        K -B+
          --+
          --+
          +++
    K
  -----+  --++
M --A--+  -C++
  -----+  --++
  ++++++  ++++"
  (let+ ((common-type (common-float-type a b))
         (m (or m (array-dimension c 0)))
         (n (or n (array-dimension c 1)))
         (k (or k (if transpose-a?
                      (array-dimension a 0)
                      (array-dimension a 1))))
         (lda (or lda (array-dimension a 1)))
         (ldb (or ldb (array-dimension b 1)))
         (ldc (or ldc (array-dimension c 1))))
    (when safe?
      (cond (transpose-a?
             (assert (<= m lda))
             (assert (<= (* lda k) (array-total-size a))))
            (t
             (assert (<= k lda))
             (assert (<= (* m lda) (array-total-size a)))))
      (cond (transpose-b?
             (assert (<= k ldb))
             (assert (<= (* ldb n) (array-total-size b))))
            (t
             (assert (<= n ldb))
             (assert (<= (* k ldb) (array-total-size b)))))
      (assert (<= n ldc))
      (assert (<= (* m ldc) (array-total-size c))))
    ;; here C=AB <=> C^T=B^T A^T, so in the argument list, A and B are
    ;; interchanged
    (blas-call ("gemm" common-type c)
      (&char (if transpose-b? #\C #\N))
      (&char (if transpose-a? #\C #\N))
      (&integers n m k) (&atom alpha)
      (&array-in b) (&integer ldb)
      (&array-in a) (&integer lda)
      (&atom beta) (&array-in/out (:input c) ())
      (&integer ldc))))

(defun scal! (alpha x &key n (incx 1) (safe? *safe?*))
  "X = alpha * X."
  (let ((type (array-float-type x))
        (n (cond (n
                  (when safe?
                    (assert (<= (* n incx) (array-total-size x))))
                  n)
                 (t (floor (array-total-size x) incx)))))
    (blas-call ("scal" type x)
      (&integer n) (&atom alpha)
      (&array-in/out (:input x) ()) (&integer incx))))

(defun axpy! (alpha x y &key n (incx 1) (incy 1) (safe? *safe?*))
  (let ((common-type (common-float-type x y))
        (n (cond (n
                  (when safe?
                    (assert (<= (* n incx) (array-total-size x)))
                    (assert (<= (* n incy) (array-total-size y))))
                  n)
                 (t (min (floor (array-total-size x) incx)
                         (floor (array-total-size y) incy))))))

    (blas-call ("axpy" common-type y)
      (&integer n) (&atom alpha)
      (&array-in/out (:input x) ()) (&integer incx)
      (&array-in/out (:input y) ()) (&integer incy))))

(defun copy! (x y &key n (incx 1) (incy 1) (safe? *safe?*))
  (let ((type (common-float-type x y))
        (n (cond (n
                  (when safe?
                    (assert (<= (* n incx) (array-total-size x)))
                    (assert (<= (* n incy) (array-total-size y))))
                  n)
                 (t (min (floor (array-total-size x) incx)
                         (floor (array-total-size y) incy))))))
    (blas-call ("copy" type y)
      (&integer n) (&array-in x) (&integer incx)
      (&array-in/out (:input y) ()) (&integer incy))))
