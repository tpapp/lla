;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:us-ascii -*-

(in-package :lla-unit-tests)
(in-readtable lla:v-syntax)


(deftestsuite basic-tests (lla-unit-tests)
  ())

;;;;
;;;; utilities
;;;;

(addtest (basic-tests)
  ;; This is basically testing that your system is ASCII.  If this
  ;; fails, then you should consider buying a more recent computer.
  ascii-characters
  (lla::with-fortran-atoms ((:char n% #\N)
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

(addtest (basic-tests)
  nv-input
  (iter 
    (for (source destination) :in (coercible-pairs-list))
    (ensure
     (test-nv-input-readonly source destination))))

;; nv-input-copied

(addtest (basic-tests)
  nv-input-copied
  (iter 
    (for (source destination) :in (coercible-pairs-list))
    (ensure
     (test-nv-input-copied source destination))))

;; nv-output

(addtest (basic-tests)
  vector-output
  (ensure (test-vector-output :integer))
  (ensure (test-vector-output :single))
  (ensure (test-vector-output :double))
  (ensure (test-vector-output :complex-single))
  (ensure (test-vector-output :complex-double)))


;;;; diagonal

(addtest (basic-tests)
  diagonal
  (ensure (== (diagonal->matrix 
               (create-diagonal '(1 2 3 4)))
              (create-matrix 4 '(1 0 0 0
                                 0 2 0 0
                                 0 0 3 0
                                 0 0 0 4))))
  (ensure (== (matrix->diagonal 
               (create-matrix 3 '(5 0 0
                                  0 7 0)))
              (create-diagonal '(5 7)))))

(addtest (basic-tests)
  diagonal-mm
  (let ((a (create-matrix 3 '(1 2 3
                              4 5 6)))
        (left-d (create-diagonal '(7 8)))
        (right-d (create-diagonal '(9 10 11))))
    (ensure (== (mm (diagonal->matrix left-d) a)
                (mm left-d a)))
    (ensure (== (mm a (diagonal->matrix right-d))
                (mm a right-d)))))


;;;;
;;;; matrix
;;;;

;;;; !!!! tests for lazy copy mechanism working correctly and copying on demand
;;;;

(addtest (basic-tests)
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

(addtest (basic-tests)
  transpose
  ;; !!! test transpose for other classes, lower-upper triangular, symmetric, etc
  ;; !!! check for type of resulting class
  (ensure (== (transpose (make-matrix-with-seq :double 3 4))
              (create-matrix 3 '(0.0d0 1.0d0 2.0d0
                                 3.0d0 4.0d0 5.0d0
                                 6.0d0 7.0d0 8.0d0
                                 9.0d0 10.0d0 11.0d0)))))

(addtest (basic-tests)
  stack
  (ensure (== (stack-vertically 
               (create-nv '(1 2 3))
               (create-matrix 3 '(4d0 5 6 7 8 9))
               (create-diagonal '(-1 -2 -3)))
              (create-matrix 3'(1 2 3
                                4 5 6
                                7 8 9
                                -1 0 0
                                0 -2 0
                                0 0 -3)))))


;;;; set-restricted

(addtest (basic-tests)
  set-restricted-and-as
  (bind ((dense #2v(1 2 3 4))
         (hermitian (as 'hermitian-matrix dense))
         (upper-triangular (as 'upper-triangular-matrix dense))
         (lower-triangular (as 'lower-triangular-matrix dense))
         (*lift-equality-test* #'equalp))
    (set-restricted hermitian)
    (ensure-same (elements hermitian) #(1 2 2 4))
    (set-restricted upper-triangular)
    (ensure-same (elements upper-triangular) #(1 0 2 4))
    (set-restricted lower-triangular)
    (ensure-same (elements lower-triangular) #(1 3 0 4))))
