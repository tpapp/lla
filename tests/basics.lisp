;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:us-ascii -*-

(in-package :lla-tests)

(deftestsuite basic-tests (lla-tests)
  ())

;;; types
;;;
;;; !! should write some basic tests for type functions, some
;;; !! conceptual tests would be nice

;;; fortran-atoms

;; (addtest (basic-tests)
;;   fortran-atoms
;;   (let ((*lift-equality-test* #'eql))
;;     (lla::with-foreign-atoms ((char% :char #\X)
;;                               (integer% :integer 1029)
;;                               (single% :single 19.0)
;;                               (double% :double 177d0)
;;                               (complex-single% :complex-single #C(5.0 3.0))
;;                               (complex-double% :complex-double #C(42d0 119d0)))
;;       (ensure-same (lla::mem-aref* char% :char) #\X)
;;       (ensure-same (lla::mem-aref* integer% :integer) 1029)
;;       (ensure-same (lla::mem-aref* single% :single) 19.0)
;;       (ensure-same (lla::mem-aref* double% :double) 177d0)
;;       (ensure-same (lla::mem-aref* complex-single% :complex-single) #C(5.0 3.0))
;;       (ensure-same (lla::mem-aref* complex-double% :complex-double) #C(42d0 119d0)))))

;; (addtest (basic-tests)
;;   ;; This is basically testing that your system is ASCII.  If this
;;   ;; fails, then you should consider buying a more recent computer.
;;   ascii-characters
;;   (lla::with-fortran-atoms ((n% :char #\N)
;;                             (c% :char #\C)
;;                             (t% :char #\T))
;;     (ensure-same (mem-ref n% :unsigned-char) 78)
;;     (ensure-same (mem-ref c% :unsigned-char) 67)
;;     (ensure-same (mem-ref t% :unsigned-char) 84)))

;; ;;; copy-elements

;; (addtest (basic-tests)
;;   copy-elements
;;   (let ((source (make-vector-with-seq 10 :single))
;;         (destination1 (lla-array 5 :single))
;;         (destination2 (lla-array 5 :double))
;;         (*lift-equality-test* #'==))
;;     (copy-elements source 2 destination1 1 3)
;;     (copy-elements source 2 destination2 2 3)
;;     (ensure-same destination1 (clo :single 0 2 3 4 0))
;;     (ensure-same destination2 (clo :double 0 0 2 3 4))))

;; (addtest (basic-tests)
;;   copy-vector
;;   (let ((source (make-vector-with-seq 3 :single))
;;         (*lift-equality-test* #'==))
;;     (ensure-same (copy-vector source) (clo :single 0 1 2))
;;     (ensure-same (copy-vector source :double) (clo :double 0 1 2))
;;     (ensure-same (copy-vector source :complex-single) (clo :complex-single 0 1 2))
;;     (ensure-same (copy-vector source :complex-double) (clo :complex-double 0 1 2))))

;; (addtest (basic-tests)
;;   copy-columns
;;   (let ((4x4-double (make-vector-with-seq 16 :double))
;;         (4x4-single (make-vector-with-seq 16 :single))
;;         (2x2 (make-vector-with-seq 4 :single))
;;         (expected-result (elements (clo :single
;;                                     0 4 8 12 :/
;;                                     1 0 2 13
;;                                     2 1 3 14
;;                                     3 7 11 15))))
;;     (copy-columns 2 2
;;                   2x2 0 2
;;                   4x4-single 5 4)
;;     (copy-columns 2 2
;;                   2x2 0 2
;;                   4x4-double 5 4)
;;     (ensure-same 4x4-single expected-result
;;                  :test #'==)
;;     (ensure-same 4x4-double (copy-vector expected-result :double)
;;                  :test #'==)))


;; ;;; pinned-vector is tested in a separate file

;; special matrixes

(addtest (basic-tests)
  special-univariate-operation
  (let ((*lift-equality-test* #'==))
    (ensure-same (e- (upper t 2)) (upper t -2))
    (ensure-same (e/ (upper t 2)) (upper t 0.5))
    (ensure-same (e+ (upper t 2 )) (upper t 2))))

(addtest (basic-tests)
  special-bivariate-operation
  (let+ ((*lift-equality-test* #'==)
         (a (dense t
              (1 2)
              (3 4)))
         (b (dense t
              (5 7)
              (11 13)))
         ((&macrolet test (kind op)
            `(ensure-same (,op (convert-matrix ,kind a)
                               (convert-matrix ,kind b))
                          (convert-matrix ,kind (,op a b)))))
         ((&macrolet tests (kind &optional (ops '(e+ e- e*)))
            `(progn
              ,@(mapcar (lambda (op) `(test ,kind ,op)) ops)))))
    (tests 'hermitian)
    (tests 'lower)
    (tests 'upper)))

(addtest (basic-tests)
  special-bivariate-to-array
  (let+ ((*lift-equality-test* #'==)
         (a (dense t
              (1 2)
              (3 4)))
         (b (dense t
              (5 7)
              (11 13)))
         ((&macrolet test (kind op)
            `(progn
               (ensure-same (,op (convert-matrix ,kind a) b) (,op a b))
               (ensure-same (,op a (convert-matrix ,kind b)) (,op a b)))))
         ((&macrolet tests (kind &optional (ops '(e+ e- e*)))
            `(progn ,@(mapcar (lambda (op) `(test ,kind ,op)) ops)))))
    (tests 'hermitian)
    (tests 'lower)
    (tests 'upper)))

(addtest (basic-tests)
  special-mean
  (let ((a (upper 'lla-double
             (2 4)
             (0 6)))
        (b (lower 'lla-double
             (2 0)
             (4 6)))
        (*lift-equality-test* #'==))
    (ensure-same (mean (list (e* 2 a) (e- a))) (e/ a 2)
                 :ignore-multiple-values? t)
    (ensure-same (mean (list (e* 2 b) (e- b))) (e/ b 2)
                 :ignore-multiple-values? t)
    (ensure-same (mean (list a b)) (mean (list (as-array a) (as-array b)))
                 :ignore-multiple-values? t)))

(addtest (basic-tests)
  wrapped-stack
  (let* ((a (diag 'double-float 1 2))
         (b (upper 'double-float
              (3 5)
              (% 7)))
         (c (vec 'double-float 11 13))
         (*lift-equality-test* #'array=))
    (ensure-same (stack nil :v a b c)
                 (dense 'double-float
                   (1 0)
                   (0 2)
                   (3 5)
                   (0 7)
                   (11 13)))
    (ensure-same (stack nil :h a b c)
                 (dense 'double-float
                   (1 0 3 5 11)
                   (0 2 0 7 13)))))
