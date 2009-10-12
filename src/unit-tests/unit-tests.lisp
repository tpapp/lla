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
  (let ((data (numeric-vector-data numeric-vector)))
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
  (make-numeric-vector length lla-type
		       (iter
			 (repeat length)
			 (collecting (random random-arg)))))

(defun test-nv-input (source-type destination-type &key
		      (report-p t) (length 50))
  "Generate a random vector copy (and convert), call nv-input and
check equality.  Tests with-nv-input (and of course nv-copy).  Types
are LLA-types.  Random integers ensure that all coercions are valid."
  (let ((nv (make-random-numeric-vector length source-type 100)))
    (with-nv-input (nv pointer destination-type)
      (numeric-vector-pointer= nv pointer destination-type :report-p
			       report-p))))
;; TESTS

;; utilities

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

;; types

;; !!! should write some basic tests for type functions, especially
;;     coercible-p and smallest-common-destination-type -- Tamas

;; nv-input

(addtest (lla-unit-tests)
  nv-input
  (for-coercible-pairs (source destination)
    (ensure
     (test-nv-input source destination))))

;; linear algebra !!! incorporate stuff below



(defparameter *a* (make-matrix 2 2 :double :initial-contents '(1 2 3 4)))
(defparameter *x* (make-matrix 2 1 :single :initial-contents '(5 6)))

(eigen *a* :check-real-p t :vectors-p t)


;;;;
;;;; output from R
;;;;
;; > a <- matrix(1:4,2,2,byrow=TRUE)
;; > a
;;      [,1] [,2]
;; [1,]    1    2
;; [2,]    3    4
;; > eigen(a)
;; $values
;; [1]  5.3722813 -0.3722813
;; $vectors
;;            [,1]       [,2]
;; [1,] -0.4159736 -0.8245648
;; [2,] -0.9093767  0.5657675

(defparameter *b* (mm *a* *x*))
(solve *a* *b*)

(solve (lu *a*) *b*)

(least-squares *a* *b*)
