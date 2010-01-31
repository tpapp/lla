;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla-unit-tests)

;;;; linear algebra

(addtest (lla-unit-tests)
  mm-solve-lu
  (bind ((a (create-matrix 2 '(1 2 3 4)))
         (x (create-nv '(5 6)))
         (b (mm a x))
         (a-lu (lu a))
         (x-solve (solve a b))
         (x-solve-lu (solve a-lu b)))
    (ensure (== b (create-nv '(17d0 39d0))))
    (ensure (== x-solve x))
    (ensure (== x-solve-lu x))))

(addtest (lla-unit-tests)
  invert
  (let ((m (create-matrix 2 '(1 2 3 4))))
    (flet ((invert (kind)
             (invert (copy-matrix m :kind kind))))
      (ensure (== (invert :dense) ; also tests (invert lu)
                  (create-matrix 2 '(-2 1 1.5 -0.5) :lla-type :double)))
      (ensure (== (invert :upper-triangular)
                  (create-matrix 2 '(1 -0.5 0 0.25) :kind :upper-triangular 
                                 :lla-type :double)))
      (ensure (== (invert :lower-triangular)
                  (create-matrix 2 '(1 0 -0.75 0.25) :kind :lower-triangular
                                 :lla-type :double))))))

(addtest (lla-unit-tests)
  eigen
  (bind ((a (create-matrix 2 '(1 2 3 4)))
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
  (bind ((a (create-matrix 2 '(1 2 3 4)))
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
  (bind ((x (create-matrix 2 '(23 23
                               22 21
                               25 20
                               29 32
                               24 29)))
         (y (create-nv '(67 63 65 94 84)))
         ((:values beta qr ss) (least-squares x y))
         (raw-var (least-squares-xxinverse qr))
         (variance (bind (((row col) (xdims x))
                          (degf (- row col))
                          (s (/ ss degf)))
                     (x* raw-var s))))
    (ensure (== (create-nv '(0.7633278 2.2350028)) beta 1d-2))
    (ensure-same ss 6.724986 :test #'approx=)
    (ensure (== variance (create-matrix 2 '(0.04035491 -0.03885797
                                            -0.03885797  0.03810950))))))

(addtest (lla-unit-tests)
  cholesky                       ; also tests hermitian factorizations
  (let* ((a (create-matrix 3 '(2 -1 0
                               -1 2 -1
                               0 -1 2)
                           :kind :hermitian))
         (r (create-matrix 3 '(1.414214 0.000000 0.0000000
                               -0.7071068 1.2247449 0.0000000
                               0.000000 -0.8164966 1.1547005)
                           :kind :upper-triangular))
         (cr (cholesky a))
         (cl (cholesky a :L))
         (b (create-nv '(5 7 13)))
         (a\b (solve a b))
         (a\1 (invert a)))
    (ensure (== (factorization-component cr :U)
                r))
    (ensure (== (factorization-component cr :L)
                (transpose r)))
    (ensure (== (factorization-component cl :U)
                r))
    (ensure (== (factorization-component cl :L)
                (transpose r)))
    (ensure (== a\b (solve cr b)))
    (ensure (== a\b (solve cl b)))
    (ensure (== a\1 (invert cr)))
    (ensure (== a\1 (invert cl)))))

(addtest (lla-unit-tests)
  svd
  (bind ((a (create-matrix 2 '(1d0 2 3 4)))
         ((:values d u vt) (svd a :all :all)))
    (ensure (== d (create-diagonal '(5.4649857 0.3659662))))
    (ensure (== u (create-matrix 2 '(-0.4045536 -0.9145143
                                     -0.9145143  0.4045536))))
    (ensure (== vt (create-matrix 2 '(-0.5760484 -0.8174156
                                      0.8174156 -0.5760484))))))

(addtest (lla-unit-tests)
  svd-rectangular
  (bind ((a (create-matrix 2 '(1 2 3 4 5 6)))
         ((:values d u vt) (svd a :singular :singular)))
    (ensure (== d (create-diagonal '(9.5255181 0.5143006))))
    (ensure (== u (create-matrix 2 '(0.2298477  0.8834610
                                     0.5247448  0.2407825
                                     0.8196419 -0.4018960))))
    (ensure (== vt (create-matrix 2 '(0.6196295 0.7848945
                                      -0.7848945 0.6196295))))))

