;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla-unit-tests)

(deftestsuite linear-algebra-tests (lla-unit-tests)
  ())

;;;; linear algebra

(addtest (linear-algebra-tests)
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

(addtest (linear-algebra-tests)
  mm-solve-triangular
  (bind ((u (clo :upper-triangular 
                 1 2 :/
                 0 3))
         (l (clo :lower-triangular
                 1 0 :/
                 2 3))
         (b (clo :dense
                 5 6 :/
                 7 8))
         (x-u (solve u b))
         (x-l (solve l b)))
    (ensure-same (mm u x-u) b :test #'==)
    (ensure-same (mm l x-l) b :test #'==)))

(addtest (linear-algebra-tests)
  invert
  (let ((m (create-matrix 2 '(1 2 3 4))))
    (flet ((invert (kind)
             (invert (copy-matrix% m :kind kind :copy-p t))))
      (ensure (== (invert :dense) ; also tests (invert lu)
                  (create-matrix 2 '(-2 1 1.5 -0.5) :lla-type :double)))
      (ensure (== (invert :upper-triangular)
                  (create-matrix 2 '(1 -0.5 0 0.25) :kind :upper-triangular 
                                 :lla-type :double)))
      (ensure (== (invert :lower-triangular)
                  (create-matrix 2 '(1 0 -0.75 0.25) :kind :lower-triangular
                                 :lla-type :double))))))

(addtest (linear-algebra-tests)
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

(addtest (linear-algebra-tests)
  eigen-single
  ;; necessary to test for single type because eigen is so hairy
  (bind ((a (create-matrix 2 '(1 2 3 4) :lla-type :single))
         ((:values eigenvalues eigenvectors)
          (eigen a :vectors-p t :check-real-p t))
         (order (xorder eigenvalues #'<)))
    ;; we order by eigenvector magnitude
    (ensure-same (slice eigenvalues order)
                 #(-0.3722813 5.3722813) :test #'x=)
    (ensure-same (slice eigenvectors :all order)
                 #2A((-0.8245648 -0.4159736)
                     (0.5657675 -0.9093767)) :test #'x=)))

(addtest (linear-algebra-tests)
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

(addtest (linear-algebra-tests)
  least-squares
  (bind ((x (create-matrix 2 '(23 23
                               22 21
                               25 20
                               29 32
                               24 29)))
         (y (create-nv '(67 63 65 94 84)))
         ((:values beta qr ss) (least-squares y x))
         (raw-var (least-squares-xx-inverse qr))
         (variance (bind (((row col) (xdims x))
                          (degf (- row col))
                          (s (/ ss degf)))
                     (x* raw-var s))))
    (ensure (== (create-nv '(0.7633278 2.2350028)) beta 1d-2))
    (ensure-same ss 6.724986 :test #'approx=)
    (ensure (== variance (create-matrix 2 '(0.04035491 -0.03885797
                                            -0.03885797  0.03810950)
                                        :kind :hermitian)))))

(addtest (linear-algebra-tests)
  constrained-least-squares
  ;;  Taken from the LAPACK documentation
  (bind ((x (clo 1 1 1 1 :/
                 1 3 1 1
                 1 -1 3 1
                 1 1 1 3
                 1 1 1 -1))
         (y (clo 2 1 6 3 1))
         (z (clo 1 1 1 -1 :/
                 1 -1 1 1
                 1 1 -1 1))
         (w (clo 1 3 -1)))
    (ensure-same (constrained-least-squares y x z w)
                 (clo 0.5 -0.5 1.5 0.5) :test #'==)))

(addtest (linear-algebra-tests)
  cholesky                       ; also tests hermitian factorizations
  (let* ((a (create-matrix 3 '(2 -1 0
                               -1 2 -1
                               0 -1 2)
                           :kind :hermitian))
         (r (create-matrix 3 '(1.414214 0.000000 0.0000000
                               -0.7071068 1.2247449 0.0000000
                               0.000000 -0.8164966 1.1547005)
                           :kind :lower-triangular))
         (cr (cholesky a))
         (cl (cholesky a :L))
         (b (create-nv '(5 7 13)))
         (a\b (solve a b))
         (a\1 (invert a)))
    (ensure (== (component cr :L)
                r))
    (ensure (== (component cr :U)
                (transpose r)))
    (ensure (== (component cl :L)
                r))
    (ensure (== (component cl :U)
                (transpose r)))
    (ensure (== a\b (solve cr b)))
    (ensure (== a\b (solve cl b)))
    (ensure (== a\1 (invert cr)))
    (ensure (== a\1 (invert cl)))))

(addtest (linear-algebra-tests)
  svd
  (bind ((a (create-matrix 2 '(1d0 2 3 4)))
         (d-true (clo :diagonal 5.4649857 0.3659662))
         (u-true (clo -0.4045536 -0.9145143 :/
                     -0.9145143  0.4045536))
         (vt-true (clo -0.5760484 -0.8174156 :/
                       0.8174156 -0.5760484))
         ((:values d u vt) (svd a :left :all :right :all))
         ((:values d2 u2 vt2) (svd a)))
    ;; all, all
    (ensure (== d d-true))
    (ensure (== u u-true))
    (ensure (== vt vt-true))
    ;; none, none
    (ensure (== d2 d-true))
    (ensure (not u2))
    (ensure (not vt2))))

(addtest (linear-algebra-tests)
  svd-rectangular
  (bind ((a (create-matrix 2 '(1 2 3 4 5 6)))
         ((:values d u vt) (svd a :left :singular :right :singular)))
    (ensure (== d (clo :diagonal 9.5255181 0.5143006)))
    (ensure (== u (clo 0.2298477  0.8834610 :/
                       0.5247448  0.2407825
                       0.8196419 -0.4018960)))
    (ensure (== vt (clo 0.6196295 0.7848945 :/
                        -0.7848945 0.6196295)))))

(addtest (linear-algebra-tests)
  tr
  (ensure-same (tr (clo 1 2 :/
                        3 4)) 5d0)
  (ensure-same (tr (clo :hermitian
                        1 2 :/
                        0 9)) 10d0)
  (ensure-same (tr (clo 1 2 :/
                        3 4)) 5d0)
  (ensure-same (tr (clo :upper-triangular
                        1 2 :/
                        3 4)) 5d0)
  (ensure-same (tr (clo :diagonal 2 15)) 17d0))

(addtest (linear-algebra-tests)
  rank
  (ensure-same (rank (clo 1 1 :/ 1 1)) 1)
  (ensure-same (rank (clo 2 4 1 3 :/
                          -1 -2 1 0
                          0 0 2 2
                          3 6 2 5)) 2))

(addtest (linear-algebra-tests)
  det
  (let ((*lift-equality-test* (x~= 1e-4)))
    (ensure-same (det (clo :upper-triangular
                           1 2 :/
                           0 4))
                 4)
    (ensure-same (det (clo :lower-triangular
                           7 0 :/
                           8 12))
                 (* 7 12))
    (ensure-same (det (mm t (clo 1 2 :/
                                 3 4)))
                 4)))
