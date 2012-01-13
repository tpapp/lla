;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla)

(defun test-mm (&optional (n 1000000))
  (let* ((a (clo :double 1 2 :/ 3 4))
        (product a))
    (dotimes (i n)
      (setf product (mm a product))
      (when (zerop (mod i 100))
        (setf product a)))))

;;; load-time-value: 12.76
;;; without: 24.79
;;; inline lookup: 22.81 
;;; foreign-funcall: 13.218

(time (test-mm))


;; Hello Tamas.  The folks here have been using LLA and recently asked me
;; to run some simple benchmarks to see how fast it is for simple things
;; like multiplying a matrix of double-floats by a vector.  I'm using a
;; 64-bit Intel workstation running Ubuntu Linux.  In the benchmarks
;; below, DOIT is the Lisp function:

(defun doit (a v)
  (declare (optimize (compilation-speed 0) (debug 0) (safety 0) (space 0)
                     (speed 3)))
  (declare (type (simple-array double-float (* *)) a)
           (type (simple-array double-float (*)) v))
  (destructuring-bind (rows columns) (array-dimensions a)
    (declare (type fixnum rows columns))
    (let ((result (make-array rows :element-type 'double-float)))
      (dotimes (row rows)
        (let ((sum 0d0))
          (dotimes (column columns)
            (incf sum (* (aref a row column) (aref v column))))
          (setf (aref result row) sum)))
      result)))

;; For my timings, I create two arrays holding random double floats and
;; check that DOIT and LLA:MMM produce the same results when multiplying
;; *A* by *V*.  Then I time the calls:

(defparameter *a* (generate-array '(50 201) (lambda () (random 1d0)) 'double-float))
(defparameter *v* (generate-array '(201 300)  (lambda () (random 1d0)) 'double-float))
(e- (lla:mmm *a* *v*) (doit *a* *v*))
(time (progn (loop repeat 10000 do (lla:mm *a* *v*)) nil))
(time (progn (loop repeat 10000 do (doit *a* *v*)) nil))

(require :sb-sprof)

(sb-sprof:with-profiling (:loop nil :report :flat)
  (time (progn (loop repeat 1000 do (lla:mm *a* *v*)) nil)))

(trace lla-to-lisp-type*)

()
(defparameter *a* (generate-array '(5 6) (lambda () (random 1d0)) 'double-float))
(defparameter *v* (generate-array 6  (lambda () (random 1d0)) 'double-float))

;;; compare to waaf-cffi

(asdf:load-system "waaf-cffi")

(defun copy-waaf (array)
  (waaf:with-array-as-foreign-pointer (array pointer :double)
    :copy-from-foreign nil)
  )

(defun copy-lla (array)
  (with-pinned-array (pointer array :double nil nil nil nil)
    ))

(time (progn (loop repeat 1000 do (copy-waaf *a*))))
(time (progn (loop repeat 1000 do (copy-lla *a*))))


(lla:mm *a* *v*)


(defparameter cl-user:*lla-configuration*
  '(:efficiency-warnings (:array-type :array-conversion)))

(let ((lla:*lla-efficiency-warning-array-type* t)
      (lla:*lla-efficiency-warning-array-conversion* t))
  (code that you want to check))

Here is a table comparing relative run times for tests like the above
using LLA:MMM (with BLAS), DOIT, and matlisp (both BLAS and ATLAS
versions):

    LLA:MMM/BLAS  17.2
    DOIT          .204
    matlisp/BLAS  .126
    matlisp/ATLAS .098


;;;; matrix copy/transpose

(defun transpose-to-memory1 (matrix pointer)
  "Specialized for double-float."
  (declare (optimize speed (safety 0)))
  (check-type matrix (simple-array double-float (* *)))
  (destructuring-bind (nrow ncol) (array-dimensions matrix)
    (let ((index 0))
      (declare (type fixnum index)
               (type (simple-array double-float (* *)) matrix))
      (loop for col-index fixnum below ncol do
        (loop for row-index fixnum below nrow do
          (setf (cffi:mem-aref pointer :double index)
                (aref matrix row-index col-index))
          (incf index))))))

(defparameter *matrix* (make-array '(100 100) :element-type 'double-float))

(time 
 (cffi:with-foreign-pointer (pointer (* (cffi:foreign-type-size :double)
                                        (array-total-size *matrix*)))
   (loop repeat 10000 do
         (transpose-to-memory1 *matrix* pointer))))

(require :sb-sprof)

(time 
 (cffi:with-foreign-pointer (pointer (* (cffi:foreign-type-size :double)
                                        (array-total-size *matrix*)))
   (loop repeat 10000 do
         (lla::transpose-matrix-to-memory *matrix* pointer lla::+double+))))


(sb-sprof:with-profiling (:loop nil :report :flat)
  (cffi:with-foreign-pointer (pointer (* (cffi:foreign-type-size :double)
                                         (array-total-size *matrix*)))
    (loop repeat 100000 do
         (lla::transpose-matrix-to-memory *matrix* pointer lla::+double+))))



