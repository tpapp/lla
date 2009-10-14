;; At the moment I am not exporting anything, so temporarily let's use
;; the LLA package.

;; (in-package :lla-unit-tests)

(in-package :lla)
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
  (let ((data (nv-data numeric-vector)))
    (dotimes (i (length data))
      (let ((nv-elt (aref data i))
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
  (make-nv length lla-type
           (iter
             (repeat length)
             (collecting (random random-arg)))))

(defun test-nv-input (source-type destination-type &key
		      (report-p t) (length 50))
  "Generate a random vector copy (and convert), use nv-input and
check equality.  Tests with-nv-input (and of course nv-copy).  Types
are LLA-types.  Random integers ensure that all coercions are valid."
  (let ((nv (make-random-numeric-vector length source-type 100)))
    (with-nv-input (nv pointer destination-type)
      (numeric-vector-pointer= nv pointer destination-type :report-p
			       report-p))))

(defun test-nv-input-copied (source-type destination-type &key
		      (report-p t) (length 50))
  "Generate a random vector copy (and convert), use nv-input-copied
and check equality.  Tests with-nv-input-copied (and of course
nv-copy).  Types are LLA-types.  Random integers ensure that all
coercions are valid."
  (let* ((nv (make-random-numeric-vector length source-type 100))
         (data-copy (copy-seq (nv-data nv))))
    (with-nv-input-copied (nv pointer destination-type)
      ;; check equality
      (and (numeric-vector-pointer= nv pointer destination-type :report-p
                                    report-p)
           ;; change location, then check that the original didn't change
           (progn
             (dotimes (i length)
               (incf (mem-aref* pointer destination-type i) 1))
             (equalp data-copy (nv-data nv)))))))

(defun test-nv-output (type &key (length 50))
  "Test with-nv-output by filling the memory are with numbers."
  (let* ((lisp-type (lla-type->lisp-type type))
         (x (with-nv-output (x length pointer type)
              (dotimes (i length)
                (setf (mem-aref* pointer type i) (coerce i lisp-type)))
              x)))
    (iter
      (for i :from 0 :below length)
      (always (= (xref x i) i)))))
           
(defparameter *allowed-difference* 1d-5
  ;; default is for catching major mistakes, use smaller value for fine points
  "Maximum allowed difference used by approx=.")

(defun approx= (a b)
  "Approximately equal, by *allowed-difference*."
  (< (abs (- a b)) *allowed-difference*))

;; TESTS

;;;;
;;;; utilities
;;;;

(addtest (lla-unit-tests)
  ;; This is basically testing that your system is ASCII.  If this
  ;; fails, then you should consider buying a more recent computer.
  with-characters
  (with-characters ((#\N n%)
		    (#\C c%)
		    (#\T t%))
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
  (for-coercible-pairs (source destination)
    (ensure
     (test-nv-input source destination))))

;; nv-input-copied

(addtest (lla-unit-tests)
  nv-input-copied
  (for-coercible-pairs (source destination)
    (ensure
     (test-nv-input-copied source destination))))

;; nv-output

(addtest (lla-unit-tests)
  nv-output
  (ensure (test-nv-output :integer))
  (ensure (test-nv-output :single))
  (ensure (test-nv-output :double))
  (ensure (test-nv-output :complex-single))
  (ensure (test-nv-output :complex-double)))

;;;;
;;;; linear algebra !!! incorporate stuff below
;;;;

(addtest (lla-unit-tests)
  mm-solve-lu
  (bind ((a (make-matrix 'dense-matrix 2 2 :initial-contents '(1 2 3 4)))
         (x (make-nv 2 :double '(5 6)))
         (b (mm a x))
         (a-lu (lu a))
         (x-solve (solve a b))
         (x-solve-lu (solve a-lu b)))
    (ensure-same b #(17d0 39d0) :test #'x=)
    (ensure-same x-solve x :test #'x=)
    (ensure-same x-solve-lu x :test #'x=)))

(addtest (lla-unit-tests)
  eigen
  (bind ((a (make-matrix 'dense-matrix 2 2 :initial-contents '(1 2 3 4)))
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
  least-squares
  (bind ((x (make-matrix 'dense-matrix 5 2 :initial-contents 
                         '(23 23 22 21 25 20 29 32 24 29)))
         (y (make-nv 5 :double '(67 63 65 94 84)))
         ((:values beta qr ss) (least-squares x y))
         (raw-var (least-squares-raw-variance qr))
         (variance (xmap 'dense-matrix 
                         (bind (((row col) (xdims x))
                                (degf (- row col))
                                (s (/ ss degf)))
                           (lambda (x)
                             (* x s)))
                         raw-var)))
    (ensure (x= #(0.7633278 2.2350028) beta))
    (ensure-same ss 6.724986 :test #'approx=)
    (ensure-same variance #2A((0.04035491 -0.03885797)
                              (-0.03885797  0.03810950))
                 :test #'x=)))

;;;; run all tests
;;;; (run-lla-tests)
