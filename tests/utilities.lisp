;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package :lla-tests)

;; SUPPORT FUNCTIONS

;; support functions

(defun random-array (type &rest dimensions)
  "Random array for testing."
  (filled-array dimensions (if (subtypep type 'complex)
                               (lambda () (coerce (complex (random 100)
                                                           (random 100))
                                                  type))
                               (lambda () (coerce (random 100) type))) type))

(defmacro with-foreign-temporary-buffer ((pointer size) &body body)
    "Allocate a buffer for SIZE complex doubles."
    `(with-foreign-pointer (,pointer (* ,size (foreign-type-size :double) 2))
       ,@body))

;;; !! following 2 functions are redundant, should be phased out
	
(defun array-at-pointer= (array pointer lla-type &key (test #'=) (report? nil))
  "Return non-nil if the elements of array are not equal (using test) to the data at
pointer, which can be of a different type, specified by lla-type.  If REPORT?,
the first discrepancy is reported to *ERROR-OUTPUT*, before returning NIL."
  (dotimes (index (array-total-size array))
    (let ((array-element (row-major-aref array index))
          (pointer-element (lla::mem-aref* pointer lla-type index)))
      (unless (funcall test array-element pointer-element)
        (when report?
          (format *error-output*
                  "at index ~a, ~a /= ~a~%" index array-element pointer-element))
        (return-from array-at-pointer= nil))))
  t)

(defun make-random-array (dimensions lla-type &optional (random-arg 100))
  "Make a vector of the given LLA type, filled with random elements.
Random is called with random-arg, the default is an integer int order
to facilitate comparing with = in case of type conversions, or if you
are using integer for lla-type."
  (aprog1 (make-array* dimensions lla-type)
    (dotimes (index (array-total-size it))
      (setf (row-major-aref it index)
            (coerce* (random random-arg) lla-type)))))

(defparameter *allowed-difference* 1d-5
  ;; default is for catching major mistakes, use smaller value for fine points
  "Maximum allowed difference used by approx=.")

;; (defun make-vector-with-seq (n lla-type)
;;   "Return a numeric-vector of type LLA-TYPE, holding the integers 0...n-1."
;;   (aprog1 (lla-array n lla-type)
;;     (dotimes (index n)
;;       (setf (aref it index) (coerce* index lla-type)))))

;; (defun make-matrix-with-seq (nrow ncol lla-type)
;;   "Return a dense matrix of type LLA-TYPE, holding the integers
;; 0...n-1 in column-major order."
;;   (lla::make-matrix% nrow ncol (make-vector-with-seq (* nrow ncol) lla-type)))

(deftestsuite utilities-tests (lla-tests)
  ())

(addtest (utilities-tests)
  approx=-test
  (let ((*==-tolerance* 1))
    (ensure (== 0 10))
    (ensure (not (== -1 10)))))

(addtest (utilities-tests)
  ==-test
  (let ((*lift-equality-test* #'==))
    (ensure-same (clo 1 :/) (clo 1 :/))
    (ensure-different (clo 1 :/) (clo 2 :/))))

;; (addtest (utilities-tests)
;;   make-vector-or-matrix-test
;;   (let ((*lift-equality-test* #'==))
;;     (ensure-same (make-vector-with-seq 3 :single) (clo :single 0 1 2))
;;     (ensure-same (make-matrix-with-seq 2 2 :complex-double)
;;                  (clo :complex-double
;;                       0 2 :/
;;                       1 3))))
