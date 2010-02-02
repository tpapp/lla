;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package :lla-unit-tests)

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
    (with-nv-input ((nv) pointer destination-type)
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
    (with-nv-input ((nv :copied) pointer destination-type)
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
