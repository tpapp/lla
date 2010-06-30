;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package :lla-unit-tests)

;; SUPPORT FUNCTIONS

(defun vector-at-pointer= (vector pointer lla-type
                           &key (test #'=) (report-p nil))
  "Return non-nil if the elements of vector are not equal \(using
test) to the data lying at pointer, which can be of a different type,
specified by lla-type.  If report-p, the first discrepancy is reported
to *error-output*, before returning nil."
  (dotimes (i (length vector))
    (let ((elt (aref vector i))
          (ptr-elt (mem-aref* pointer lla-type i)))
      (unless (funcall test elt ptr-elt)
        (when report-p
          (format *error-output*
                  "at index ~a, ~a /= ~a~%" i elt ptr-elt))
        (return-from vector-at-pointer= nil))))
  t)
	
(defun make-random-vector (length lla-type &optional
                           (random-arg 100))
  "Make a vector of the given LLA type, filled with random elements.
Random is called with random-arg, the default is an integer int order
to facilitate comparing with = in case of type conversions, or if you
are using integer for LLA-type."
  (aprog1 (lla-vector length lla-type)
    (dotimes (i length)
      (setf (aref it i) (coerce* (random random-arg) lla-type)))))


(defparameter *allowed-difference* 1d-5
  ;; default is for catching major mistakes, use smaller value for fine points
  "Maximum allowed difference used by approx=.")

(defun approx= (a b)
  "Approximately equal, by *allowed-difference*."
  (< (abs (- a b)) *allowed-difference*))

(defun == (a b &optional (eps *allowed-difference*))
  "`Strict' equality: A and B have the same class, same dimensions,
  elements differ at most by EPS.  Error is signalled on mismatch."
  (labels ((elementwise= (a-elements b-elements)
             (assert (= (length a-elements) (length b-elements)) ()
                     "elements don't have equal length")
             (assert (equalp (array-element-type a-elements)
                             (array-element-type b-elements)) ()
                             "elemens don't have the same type")
             (iter
               (for a :in-vector a-elements)
               (for b :in-vector b-elements)
               (let ((diff (abs (- a b))))
                 (assert (< diff eps) () "elements differ by at least ~A." diff)))
             t)
           (elements= ()
             (when (typep a 'dense-matrix-like)
               (set-restricted a))
             (when (typep b 'dense-matrix-like)
               (set-restricted b))
             (if (typep a 'lla::elements%)
                 (elementwise= (elements a) (elements b))
                 (elementwise= a b))))
    (assert (equal (class-of a) (class-of b)) ()
            "objects don't have the same class")
    (when (typep a 'dense-matrix)
      (assert (and (= (nrow a) (nrow b)) (= (ncol a) (ncol b))) ()
              "the matrices don't have the same dimension"))
    (elements=)))

(defun make-vector-with-seq (n lla-type)
  "Return a numeric-vector of type LLA-TYPE, holding the integers 0...n-1."
  (aprog1 (lla-vector n lla-type)
    (dotimes (index n)
      (setf (aref it index) (coerce* index lla-type)))))

(defun make-matrix-with-seq (nrow ncol lla-type)
  "Return a dense matrix of type LLA-TYPE, holding the integers
0...n-1 in column-major order."
  (lla::make-matrix% nrow ncol (make-vector-with-seq (* nrow ncol) lla-type)))
