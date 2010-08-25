;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:us-ascii -*-

(in-package :lla-unit-tests)

(deftestsuite basic-tests (lla-unit-tests)
  ())

(addtest (basic-tests)
  ;; This is basically testing that your system is ASCII.  If this
  ;; fails, then you should consider buying a more recent computer.
  ascii-characters
  (lla::with-fortran-atoms ((n% :char #\N)
                            (c% :char #\C)
                            (t% :char #\T))
    (ensure-same (mem-ref n% :unsigned-char) 78)
    (ensure-same (mem-ref c% :unsigned-char) 67)
    (ensure-same (mem-ref t% :unsigned-char) 84)))

;;; utilities
;;;
;;; !! nothing tested at the moment, all functions are very basic


;;; types
;;;
;;; !! should write some basic tests for type functions, some
;;; !! conceptual tests would be nice


;;; fortran-atoms

(addtest (basic-tests)
  fortran-atoms
  (let ((*lift-equality-test* #'eql))
    (with-fortran-atoms ((char% :char #\X)
                         (integer% :integer 1029)
                         (single% :single 19.0)
                         (double% :double 177d0)
                         (complex-single% :complex-single #C(5.0 3.0))
                         (complex-double% :complex-double #C(42d0 119d0)))
      (ensure-same (lla::mem-aref* char% :char) #\X)
      (ensure-same (lla::mem-aref* integer% :integer) 1029)
      (ensure-same (lla::mem-aref* single% :single) 19.0)
      (ensure-same (lla::mem-aref* double% :double) 177d0)
      (ensure-same (lla::mem-aref* complex-single% :complex-single) #C(5.0 3.0))
      (ensure-same (lla::mem-aref* complex-double% :complex-double) #C(42d0 119d0)))))


;;; copy-elements

(addtest (basic-tests)
  copy-elements
  (let ((source (make-vector-with-seq 10 :single))
        (destination1 (lla-array 5 :single))
        (destination2 (lla-array 5 :double))
        (*lift-equality-test* #'==))
    (copy-elements source 2 destination1 1 3)
    (copy-elements source 2 destination2 2 3)
    (ensure-same destination1 (clo :single 0 2 3 4 0))
    (ensure-same destination2 (clo :double 0 0 2 3 4))))

(addtest (basic-tests)
  copy-vector
  (let ((source (make-vector-with-seq 3 :single))
        (*lift-equality-test* #'==))
    (ensure-same (copy-vector source) (clo :single 0 1 2))
    (ensure-same (copy-vector source :double) (clo :double 0 1 2))
    (ensure-same (copy-vector source :complex-single) (clo :complex-single 0 1 2))
    (ensure-same (copy-vector source :complex-double) (clo :complex-double 0 1 2))))

(addtest (basic-tests)
  copy-columns
  (let ((4x4-double (make-vector-with-seq 16 :double))
        (4x4-single (make-vector-with-seq 16 :single))
        (2x2 (make-vector-with-seq 4 :single))
        (expected-result (elements (clo :single
                                    0 4 8 12 :/
                                    1 0 2 13
                                    2 1 3 14
                                    3 7 11 15))))
    (copy-columns 2 2
                  2x2 0 2
                  4x4-single 5 4)
    (copy-columns 2 2
                  2x2 0 2
                  4x4-double 5 4)
    (ensure-same 4x4-single expected-result 
                 :test #'==)
    (ensure-same 4x4-double (copy-vector expected-result :double) 
                 :test #'==)))


;;; pinned-vector is tested in a separate file

;;; matrix-classes

(addtest (basic-tests)
  matrix-classes
  (let ((m (make-matrix 2 2 :double)))
   (setf (mref m 0 0) 6d0)
   (incf (mref m 0 1) 7d0)
   (decf (mref m 1 0) 6d0)
   (setf (mref m 1 1) 5d0)
   (ensure-same m (clo :double
                       6 7 :/
                       -6 5)
                :test #'==)
   (let ((l (copy-matrix m :kind :lower))
         (u (copy-matrix m :kind :upper))
         (h (copy-matrix m :kind :hermitian)))
     ;; we didn't call set-restricted, so we are testing mref here
     (ensure-same (mref l 0 1) 0d0)
     (ensure-same (mref u 1 0) 0d0)
     (ensure-same (mref h 1 0) 7d0)
     (let ((*lift-equality-test* #'==))
       ;; set-restricted is called by ==
       (ensure-same l (clo :lower :double
                           6 0 :/
                           -6 5))
       (ensure-same u (clo :upper :double
                           6 7 :/
                           0 5))
       (ensure-same h (clo :hermitian :double
                           6 7 :/
                           7 5))))))

;;; diagonal

(addtest (basic-tests)
  diagonal
  (let ((d1 (clo :diagonal :double 1 2 3))
        (d2 (lla::make-diagonal% (clo :double 1 2 3))))
    (ensure-same d1 d2 :test #'==)))

;;; bind-extensions

(addtest (basic-tests)
  bind-lla-vector
  (bind (((:lla-vector a :length length) (lla-array 5 :double)))
    (dotimes (i length)
      (setf (a i) (coerce* i :double)))
    (dotimes (i length)
      (incf (a i)))
    (ensure-same a (clo :double 1 2 3 4 5) :test #'==)))

(addtest (basic-tests)
  bind-lla-matrix
  (bind (((:lla-matrix a :nrow a-nrow :ncol a-ncol) (make-matrix 3 4 :double)))
    (dotimes (row a-nrow)
      (dotimes (col a-ncol)
        (setf (a (a-index row col)) (+ (* 10d0 row) col))))
    (ensure-same a (clo :double
                        0 1 2 3 :/
                        10 11 12 13
                        20 21 22 23)
                 :test #'==)))

;;; clo

(addtest (basic-tests)
  ;; only test corner cases of CLO, since pretty much everything is
  ;; using it anyway
  (ensure-error (eval '(clo :invalid-keyword 1 2 3)))
  (ensure-error (eval '(clo 1 2 :/ 3))))
