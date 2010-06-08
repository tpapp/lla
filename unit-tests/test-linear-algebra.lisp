;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla-unit-tests)

(deftestsuite linear-algebra-tests (lla-unit-tests)
  ()
  (:equality-test #'==))

;;;; linear algebra

(addtest (linear-algebra-tests)
  diagonal-mm
  (let ((a (clo 1 2 3 :/
                4 5 6))
        (left-d (clo :diagonal 7 8))
        (right-d (clo :diagonal 9 10 11)))
    (ensure-same (mm (as-matrix left-d) a)
                 (mm left-d a))
    (ensure-same (mm a (as-matrix right-d))
                 (mm a right-d))))

(addtest (linear-algebra-tests)
  mm-solve-lu
  (bind ((a (clo 1 2 :/
                 3 4))
         (x (clo :double 5 6))
         (b (mm a x))
         (a-lu (lu a))
         (x-solve (solve a b))
         (x-solve-lu (solve a-lu b))
         )
    (ensure-same b (clo :double 17d0 39d0))
    (ensure-same x-solve x)
    (ensure-same x-solve-lu x)))

(addtest (linear-algebra-tests)
  mm-solve-triangular
  (bind ((u (clo :upper
                 1 2 :/
                 0 3))
         (l (clo :lower
                 1 0 :/
                 2 3))
         (b (clo :double
                 5 6 :/
                 7 8))
         (x-u (solve u b))
         (x-l (solve l b)))
    (ensure-same (mm u x-u) b)
    (ensure-same (mm l x-l) b)))

(addtest (linear-algebra-tests)
  invert
  (let ((m (clo :double
                1 2 :/
                3 4)))
    (flet ((invert (kind)
             (invert (copy-matrix m :kind kind :copy? t))))
      (ensure-same (invert :dense)           ; also tests (invert lu)
                   (clo :double
                        -2 1 :/ 
                        1.5 -0.5))
      (ensure-same (invert :upper)
                   (clo :double :upper 
                        1 -0.5 :/
                        0 0.25))
      (ensure-same (invert :lower)
                   (clo :double :lower
                        1 0 :/
                        -0.75 0.25)))))

(addtest (linear-algebra-tests)
  eigen
  (bind ((a (clo :double
                 1 2 :/
                 3 4))
         ((:values eigenvalues eigenvectors)
          (eigen a :vectors? t :check-real? t))
         ((:values ev order) (sort-order eigenvalues #'<)))
    ;; we order by eigenvector magnitude
    (ensure-same eigenvalues
                 (clo :double -0.3722813 5.3722813))
    (ensure-same (sub eigenvectors t order)
                 (clo :double
                      -0.8245648 -0.4159736 :/
                      0.5657675 -0.9093767))))

(addtest (linear-algebra-tests)
  eigen-single
  ;; necessary to test for single type because eigen is so hairy
  (bind ((a (clo :single
                 1 2 :/
                 3 4))
         ((:values eigenvalues eigenvectors)
          (eigen a :vectors? t :check-real? nil))
         ((:values ev order) (sort-order eigenvalues #'<
                                         :key #'abs)))
    ;; we order by eigenvector magnitude
    (ensure-same eigenvalues
                 (clo :complex-single -0.3722813 5.3722813))
    (ensure-same (sub eigenvectors t order)
                 (clo :complex-single
                      -0.8245648 -0.4159736 :/
                      0.5657675 -0.9093767))))

(addtest (linear-algebra-tests)
  eigen-complex
  (bind ((a (clo :complex-double
                 #C(1 2) #C(3 4) :/
                 #C(5 6) #C(7 8)))
         ((:values eigenvalues eigenvectors)
          (eigen a :vectors? t :check-real? nil))
         ((:values ev order) (sort-order eigenvalues #'<
                                         :key #'abs)))
    ;; we order by eigenvector magnitude
    (ensure-same eigenvalues
                 (clo :complex-double
                      #C(-0.8845984 -0.7323034)
                      #C(8.8845984 10.7323034)))
    (ensure-same (sub eigenvectors t order)
                 (clo :complex-double
                      #C(0.8331344 0.0000000) #C(0.3895109 0.0355148) :/
                      #C( -0.5526350 -0.0219453) #C(0.9203369 0.0000000)))))

(addtest (linear-algebra-tests)
  hermitian
  (bind ((a (clo :double
                 1 2 :/
                 3 4))
         (aa (mm t a))
         ((:values eigenvalues eigenvectors)
          (eigen aa :vectors? t))
         ((:values ev order) (sort-order eigenvalues #'<)))
    ;; we order by eigenvector magnitude
    (ensure-same ev
                 (clo :double 0.1339313 29.8660687))
    (ensure-same (sub eigenvectors t order)
                 (clo :double
                      -0.8174156 0.5760484 :/
                      0.5760484 0.8174156))))

(addtest (linear-algebra-tests)
  least-squares
  (bind ((x (clo :double 
                 23 23 :/
                 22 21
                 25 20
                 29 32
                 24 29))
         (y (clo :double 67 63 65 94 84))
         (beta (clo :double 0.7633278 2.2350028))
         (ss 6.724986d0)
         ((:values beta1 ss1 nu1 other1) (least-squares y x :method :svd-d))
         ((:values beta2 ss2 nu2 other2) (least-squares y x :method :qr))
         (raw-var (qr-xx-inverse (getf other2 :qr)))
         (variance (e* raw-var (/ ss2 nu2)))
         (*allowed-difference* 1e-5))
    (ensure-same beta1 beta)
    (ensure-same beta2 beta)
    (ensure-same ss1 ss :test #'approx=)
    (ensure-same ss2 ss :test #'approx=)
    (ensure-same nu1 3 :test #'=)
    (ensure-same nu2 3 :test #'=)
    (ensure-same variance (clo :double :hermitian
                               0.04035491 -0.03885797 :/
                               *  0.03810950))))

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
         (w (clo 1 3 -1))
         (*allowed-difference* 1e-5))
    (ensure-same (constrained-least-squares y x z w)
                 (clo :single 0.5 -0.5 1.5 0.5))))

(addtest (linear-algebra-tests)
  cholesky                       ; also tests hermitian factorizations
  (let* ((a (clo :double :hermitian
                 2 -1 0 :/
                 -1 2 -1
                 0 -1 2))
         (r (clo :double :lower
                 1.414214 0.000000 0.0000000 :/
                 -0.7071068 1.2247449 0.0000000
                 0.000000 -0.8164966 1.1547005))
         (cr (cholesky a))
         (cl (cholesky a :L))
         (b (clo :double 5 7 13))
         (a\b (solve a b))
         (a\1 (invert a)))
    (ensure-same (component cr :L) r)
    (ensure-same (component cr :U) (transpose r))
    (ensure-same (component cl :L) r)
    (ensure-same (component cl :U) (transpose r))
    (ensure-same a\b (solve cr b))
    (ensure-same a\b (solve cl b))
    (ensure-same a\1 (invert cr))
    (ensure-same a\1 (invert cl))))

(addtest (linear-algebra-tests)
  svd
  (bind ((a (clo :double
                 1 2 :/
                 3 4))
         (d-true (clo :diagonal :double 5.4649857 0.3659662))
         (u-true (clo :double 
                      -0.4045536 -0.9145143 :/
                     -0.9145143  0.4045536))
         (vt-true (clo :double
                       -0.5760484 -0.8174156 :/
                       0.8174156 -0.5760484))
         ((:values d u vt) (svd a :left :all :right :all))
         ((:values d2 u2 vt2) (svd a)))
    ;; all, all
    (ensure-same d d-true)
    (ensure-same u u-true)
    (ensure-same vt vt-true)
    ;; none, none
    (ensure-same d2 d-true)
    (ensure (not u2))
    (ensure (not vt2))))

(addtest (linear-algebra-tests)
  svd-rectangular
  (bind ((a (clo :double 
                 1 2 :/
                 3 4
                 5 6))
         ((:values d u vt) (svd a :left :singular :right :singular)))
    (ensure-same d (clo :diagonal :double 9.5255181 0.5143006))
    (ensure-same u (clo :double
                        0.2298477  0.8834610 :/
                        0.5247448  0.2407825
                        0.8196419 -0.4018960))
    (ensure-same vt (clo :double
                         0.6196295 0.7848945 :/
                         -0.7848945 0.6196295))))

(addtest (linear-algebra-tests)
  tr
  (let ((*lift-equality-test* #'=))
   (ensure-same (tr (clo :double
                         1 2 :/
                         3 4)) 5d0)
    (ensure-same (tr (clo :hermitian
                          :double
                          1 2 :/
                          0 9)) 10d0)
    (ensure-same (tr (clo :double
                          1 2 :/
                          3 4)) 5d0)
    (ensure-same (tr (clo :upper
                          :double
                          1 2 :/
                          3 4)) 5d0)
    (ensure-same (tr (clo :double :diagonal 2 15)) 17d0)))

(addtest (linear-algebra-tests)
  rank
  (let ((*lift-equality-test* #'=))
    (ensure-same (rank (clo 1 1 :/ 1 1)) 1)
    (ensure-same (rank (clo 2 4 1 3 :/
                            -1 -2 1 0
                            0 0 2 2
                            3 6 2 5)) 2)))

(addtest (linear-algebra-tests)
  det
  (let ((*allowed-difference* 1e-4)
        (*lift-equality-test* #'approx=))
    (ensure-same (det (clo :upper
                           1 2 :/
                           0 4))
                 4)
    (ensure-same (det (clo :lower
                           7 0 :/
                           8 12))
                 (* 7 12))
    (ensure-same (det (mm t (clo 1 2 :/
                                 3 4)))
                 4)))
