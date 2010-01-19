;; At the moment I am not exporting anything, so temporarily let's use
;; the LLA package.

;; (in-package :lla-unit-tests)

(in-package :lla)
(in-readtable lla-readtable)
(use-package :lift)

;; TEST SUITES

(deftestsuite lla-unit-tests () ())

;; EXTERNAL

(defun run-lla-tests ()
  "Run all the tests for LLA."
  (run-tests :suite 'lla-unit-tests))

;; SUPPORT FUNCTIONS

(defun numeric-vector-pointer= (numeric-vector pointer lla-type
				&key (test #'=) (report-p nil))
  "Return non-nil if the elements of numeric-vector are not equal
\(using test) to the data lying at pointer, which can be of a
different type, specified by lla-type.  If report-p, the first
discrepancy is reported to *error-output*, before returning nil."
  (let ((elements (elements numeric-vector)))
    (dotimes (i (length elements))
      (let ((nv-elt (aref elements i))
	    (ptr-elt (mem-aref* pointer lla-type i)))
	(unless (funcall test nv-elt ptr-elt)
	  (when report-p
	    (format *error-output*
		    "at index ~a, ~a /= ~a~%" i nv-elt ptr-elt))
	  (return-from numeric-vector-pointer= nil))))
    t))
	
(defun make-random-numeric-vector (length lla-type &optional
				   (random-arg 100))
  "Make a numeric-vector of the given type, filled with random
elements.  Random is called with random-arg, the default is an integer
int order to facilitate comparing with = in case of type conversions,
or if you are using integer for LLA-type."
  (let ((nv (make-nv length lla-type))
        (lisp-type (lla-type->lisp-type lla-type)))
    (dotimes (i length)
      (setf (xref nv i) (coerce (random random-arg) lisp-type)))
    nv))

(defun test-nv-input-readonly (source-type destination-type &key
                               (report-p t) (length 50))
  "Generate a random vector, copy (and convert), use
WITH-NV-INPUT-READONLY and check equality.  Types are LLA-types.
Random numbers are integers, which ensures that all coercions are
valid."
  (let ((nv (make-random-numeric-vector length source-type 100)))
    (with-nv-input-only (nv pointer destination-type nil)
      (numeric-vector-pointer= nv pointer destination-type :report-p
			       report-p))))

(defun test-nv-input-copied (source-type destination-type &key
                             (report-p t) (length 50))
  "Generate a random vector, copy (and convert), use
WITH-NV-INPUT-COPIED and check equality.  Tests WITH-NV-INPUT-COPIED
and COPY-ELEMENTS.  Types are LLA-types.  Random numbers are integers,
which ensures that all coercions are valid."
  (let* ((nv (make-random-numeric-vector length source-type 100))
         (elements-copy (copy-seq (elements nv))))
    (with-nv-input-only (nv pointer destination-type t)
      ;; check equality
      (and (numeric-vector-pointer= nv pointer destination-type :report-p
                                    report-p)
           ;; change location, then check that the original didn't change
           (progn
             (dotimes (i length)
               (incf (mem-aref* pointer destination-type i) 1))
             (equalp elements-copy (elements nv)))))))

(defun test-vector-output (type &key (length 50))
  "Test WITH-VECTOR-OUTPUT by filling the memory are with numbers."
  (let* ((lisp-type (lla-type->lisp-type type))
         (x (with-vector-output (x length pointer type)
              (dotimes (i length)
                (setf (mem-aref* pointer type i) (coerce i lisp-type)))
              x)))
    (iter
      (for i :from 0 :below length)
      (always (= (aref x i) i)))))
           
(defparameter *allowed-difference* 1d-5
  ;; default is for catching major mistakes, use smaller value for fine points
  "Maximum allowed difference used by approx=.")

(defun approx= (a b)
  "Approximately equal, by *allowed-difference*."
  (< (abs (- a b)) *allowed-difference*))

(defun x~= (eps)
  "Return a comparison function for xref'able elements, with Lsup norm eps."
  (lambda (a b)
    (x= a b eps)))

(defun == (a b &optional (eps *allowed-difference*))
  "`Strict' equality: A and B have the same class, same dimensions,
  elements differ at most by EPS.  Error is signalled on mismatch."
  (flet ((elements= ()
           (when (typep a 'restricted-elements)
             (set-restricted a))
           (when (typep b 'restricted-elements)
             (set-restricted b))
           (let ((a-elements (elements a))
                 (b-elements (elements a)))
             (assert (= (length a-elements) (length b-elements)) ()
                     "elements don't have equal length")
             (iter
               (for a :in-vector a-elements)
               (for b :in-vector b-elements)
               (let ((diff (abs (- a b))))
                 (assert (< diff eps) () "elements differ by at least ~A." diff )))
             t)))
    (assert (equal (class-of a) (class-of b)) ()
            "objects don't have the same class")
    (when (typep a 'dense-matrix-like)
      (assert (and (= (nrow a) (nrow b)) (= (ncol a) (ncol b))) ()
              "the matrices don't have the same dimension"))
    (elements=)))

(defun make-nv-with-seq (lla-type n)
  "Return a numeric-vector of type LLA-TYPE, holding the integers 0...n-1."
  (bind ((nv (make-nv n lla-type))
         (elements (elements nv))
         (lisp-type (xelttype nv)))
    (iter
      (for i :from 0 :below n)
      (setf (aref elements i) (coerce i lisp-type)))
    nv))

(defun make-matrix-with-seq (lla-type nrow ncol)
  "Return a dense matrix of type LLA-TYPE, holding the integers
0...n-1 in column-major order."
  (vector->matrix (make-nv-with-seq lla-type (* nrow ncol)) nrow ncol))

;; TESTS

;;;;
;;;; utilities
;;;;

(addtest (lla-unit-tests)
  ;; This is basically testing that your system is ASCII.  If this
  ;; fails, then you should consider buying a more recent computer.
  ascii-characters
  (with-fortran-atoms ((:char n% #\N)
                       (:char c% #\C)
                       (:char t% #\T))
    (ensure-same (mem-ref n% :unsigned-char) 78)
    (ensure-same (mem-ref c% :unsigned-char) 67)
    (ensure-same (mem-ref t% :unsigned-char) 84)))
;;;;
;;;; types
;;;;

;; !!! should write some basic tests for type functions, especially
;;     coercible-p and smallest-common-destination-type -- Tamas


;;;;
;;;; numeric-vector
;;;;

;; nv-input

(addtest (lla-unit-tests)
  nv-input
  (iter 
    (for (source destination) :in (coercible-pairs-list))
    (ensure
     (test-nv-input-readonly source destination))))

;; nv-input-copied

(addtest (lla-unit-tests)
  nv-input-copied
  (iter 
    (for (source destination) :in (coercible-pairs-list))
    (ensure
     (test-nv-input-copied source destination))))

;; nv-output

(addtest (lla-unit-tests)
  vector-output
  (ensure (test-vector-output :integer))
  (ensure (test-vector-output :single))
  (ensure (test-vector-output :double))
  (ensure (test-vector-output :complex-single))
  (ensure (test-vector-output :complex-double)))


;;;; diagonal

(addtest (lla-unit-tests)
  diagonal
  (ensure (== (diagonal->matrix #v:diagonal(1 2 3 4))
              #4vd(1 0 0 0
                   0 2 0 0
                   0 0 3 0
                   0 0 0 4)))
  (ensure (== (matrix->diagonal #3vd(5 0 0
                                     0 7 0))
              #v:diagonal(5 7))))

(addtest (lla-unit-tests)
  diagonal-mm
  (let ((a #3v(1 2 3
               4 5 6))
        (left-d #v:diagonal(7 8))
        (right-d #v:diagonal(9 10 11)))
    (ensure (== (mm (diagonal->matrix left-d) a)
                (mm left-d a)))
    (ensure (== (mm a (diagonal->matrix right-d))
                (mm a right-d)))))


;;;;
;;;; matrix
;;;;

;;;; !!!! tests for lazy copy mechanism working correctly and copying on demand
;;;;

(addtest (lla-unit-tests)
  create-matrix
  ;; Note: here we compare to Lisp arrays.  If we compared to LLA
  ;; objects created with the #v read macro, we could never detect
  ;; bugs in CREATE-MATRIX as it is used by the readmacro itself.
  (flet ((make (kind)
           (create-matrix 2 '(1 2 3 4) :kind kind)))
    (ensure-same (make :dense)
                 #2A((1 2) (3 4))
                 :test #'x=)
    (ensure-same (make :upper-triangular)
                 #2A((1 2) (0 4))
                 :test #'x=)
    (ensure-same (make :lower-triangular)
                 #2A((1 0) (3 4))
                 :test #'x=)
    (ensure-same (make :hermitian)
                 #2A((1 2) (2 4))
                 :test #'x=)))

(addtest (lla-unit-tests)
  transpose
  ;; !!! test transpose for other classes, lower-upper triangular, symmetric, etc
  ;; !!! check for type of resulting class
  (ensure (== (transpose (make-matrix-with-seq :double 3 4))
              #3vd:dense(0.0d0 1.0d0 2.0d0
                         3.0d0 4.0d0 5.0d0
                         6.0d0 7.0d0 8.0d0
                         9.0d0 10.0d0 11.0d0))))

(addtest (lla-unit-tests)
  stack
  (ensure (== (stack-vertically #v(1 2 3) #3v(4d0 5 6 7 8 9) #v(-1 -2 -3))
              #3v(1 2 3
                  4 5 6
                  7 8 9
                  -1 -2 -3))))


;;;;
;;;; linear algebra !!! incorporate stuff below
;;;;

(addtest (lla-unit-tests)
  mm-solve-lu
  (bind ((a #2v(1 2 3 4))
         (x #v(5 6))
         (b (mm a x))
         (a-lu (lu a))
         (x-solve (solve a b))
         (x-solve-lu (solve a-lu b)))
    (ensure (== b #v(17d0 39d0)))
    (ensure (== x-solve x))
    (ensure (== x-solve-lu x))))

(addtest (lla-unit-tests)
  invert
  (let ((m #2vs(1 2 3 4)))
    (flet ((invert (kind)
             (invert (copy-matrix m :kind kind))))
      (ensure (== (invert :dense) ; also tests (invert lu)
                  #2vs(-2 1 1.5 -0.5)))
      (ensure (== (invert :upper-triangular)
                  #2vs:upper-triangular(1 -0.5 0 0.25)))
      (ensure (== (invert :lower-triangular)
                  #2vs:lower-triangular(1 0 -0.75 0.25))))))

(addtest (lla-unit-tests)
  eigen
  (bind ((a #2vd(1 2 3 4))
         ((:values eigenvalues eigenvectors)
          (eigen a :vectors-p t :check-real-p t))
         (order (xorder eigenvalues #'<)))
    ;; we order by eigenvector magnitude
    (ensure-same (slice eigenvalues order)
                 #(-0.3722813 5.3722813) :test #'x=)
    (ensure-same (slice eigenvectors :all order)
                 #2A((-0.8245648 -0.4159736)
                     (0.5657675 -0.9093767)) :test #'x=)))

(addtest (lla-unit-tests)
  hermitian
  (bind ((a #2v(1 2 3 4))
         (aa (mm t a))
         ((:values eigenvalues eigenvectors)
          (eigen aa :vectors-p t))
         (order (xorder eigenvalues #'<)))
    ;; we order by eigenvector magnitude
    (ensure-same (slice eigenvalues order)
                 #(0.1339313 29.8660687) :test #'x=)
    (ensure-same (slice eigenvectors :all order)
                 #2A(( -0.8174156 0.5760484)
                     (0.5760484 0.8174156)) :test #'x=)))

(addtest (lla-unit-tests)
  least-squares
  (bind ((x #2v(23 23
                22 21
                25 20
                29 32
                24 29))
         (y #v(67 63 65 94 84))
         ((:values beta qr ss) (least-squares x y))
         (raw-var (least-squares-xxinverse qr))
         (variance (bind (((row col) (xdims x))
                          (degf (- row col))
                          (s (/ ss degf)))
                     (x* raw-var s))))
    (ensure (== #vd(0.7633278 2.2350028) beta 1d-2))
    (ensure-same ss 6.724986 :test #'approx=)
    (ensure (== variance #2vd(0.04035491 -0.03885797
                             -0.03885797  0.03810950)))))

(addtest (lla-unit-tests)
  cholesky
  (let* ((a #3vd:hermitian(2 -1 0
                           -1 2 -1
                            0 -1 2))
         (c (cholesky a)))
    (ensure (== (transpose (factorization-component c :R))
                 #3vd:lower-triangular(1.414214 0.000000 0.0000000
                                       -0.7071068 1.2247449 0.0000000
                                        0.000000 -0.8164966 1.1547005)))
    (ensure (== (reconstruct c) a))))

(addtest (lla-unit-tests)
  svd
  (bind ((a #2v(1d0 2 3 4))
         ((:values d u vt) (svd a :all :all)))
    (ensure (== d #vd(5.4649857 0.3659662)))
    (ensure (== u #2vd(-0.4045536 -0.9145143
                       -0.9145143  0.4045536)))
    (ensure (== vt #2vd(-0.5760484 -0.8174156
                        0.8174156 -0.5760484)))))

(addtest (lla-unit-tests)
  svd-rectangular
  (bind ((a #2v(1d0 2 3 4 5 6))
         ((:values d u vt) (svd a :singular :singular)))
    (ensure (== d #vd(9.5255181 0.5143006)))
    (ensure (== u #2vd(0.2298477  0.8834610
                       0.5247448  0.2407825
                       0.8196419 -0.4018960)))
    (ensure (== vt #2vd(0.6196295 0.7848945
                        -0.7848945 0.6196295)))))


;;;; run all tests
;;;; (run-lla-tests)
