;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla-tests)

(defsuite linear-algebra-suite (tests))

;;;; linear algebra

(deftest mm (linear-algebra-suite)
  (let ((a (mx 'lla-double
             (1 2)
             (3 4)
             (5 6)))
        (b2 (vec 'lla-double 1 2))
        (b3 (vec 'lla-double 1 2 3))
        (i2 (mx 'lla-double
              (1 0)
              (0 1)))
        (i3 (mx 'lla-double
              (1 0 0)
              (0 1 0)
              (0 0 1))))
    (assert-equality #'num= (mm a i2) a)
    (assert-equality #'num= (mm i3 a) a)
    (assert-equality #'num= (mm a b2) (vec 'lla-double 5 11 17))
    (assert-equality #'num= (mm b3 a) (vec 'lla-double 22 28))
    (assert-condition error (mm a b3))
    ;; (assert-condition error (mm b3 b3))
    (assert-condition error (mm b2 b3))))

(deftest mm-dot (linear-algebra-suite)
  (let+ ((a (vec 'lla-double 2 3 5))
         (b (vec 'lla-complex-double 1 #C(2 1) 3))
         ((&flet mm-vec (a b)
            (sum (map 'vector (lambda (a b)
                                (* a (conjugate b)))
                      a b))))
         ((&flet test (a1 b1 a2 b2)
            (assert-equality #'num= (mm a1 b1)
                (mm-vec a2 b2)))))
    (test a b a b)
    (test b a b a)
    (test a t a a)
    (test b t b b)
    (test t a a a)
    (test t b b b)))

(deftest outer (linear-algebra-suite)
  (let+ ((a (vec 'lla-double 2 3))
         (b (vec 'lla-complex-double 1 #C(2 1) 9))
         ((&flet outer2 (a b &optional hermitian?)
            (let ((result (mm (aops:reshape-col a)
                              (aops:reshape-row (aops:each #'conjugate b)))))
              (if hermitian?
                  (hermitian-matrix result)
                  result)))))
    (assert-equality #'num= (outer a b) (outer2 a b))
    (assert-equality #'num= (outer a t) (outer2 a a t))
    (assert-equality #'num= (outer t a) (outer2 a a t))
    (assert-equality #'num= (outer (vec 'lla-complex-double 1 #C(2 1) 9) t)
        (hermitian-mx 'lla-complex-double
          (1)
          (#C(2 1) 5)
          (9 #C (18 -9) 81)))))

(deftest mm-diagonal (linear-algebra-suite)
  (let ((a (mx 'lla-double
             (1 2)
             (3 4)
             (5 6)))
        (d2 (diagonal-mx 'lla-double 2 3))
        (d3 (diagonal-mx 'lla-double 2 3 5))
        (d3*a (mx 'lla-double
                (2 4)
                (9 12)
                (25 30)))
        (a*d2 (mx 'lla-double
                (2 6)
                (6 12)
                (10 18))))
    (assert-equality #'num= (mm a d2) a*d2)
    (assert-equality #'num= (mm d3 a) d3*a)
    (assert-condition error (mm a d3))
    (assert-condition error (mm d2 a))
    (assert-equality #'num= (mm d3 d3) (e* d3 d3))
    (assert-equality #'num= (mm d3 t) (e* d3 d3))
    (assert-equality #'num= (mm t d3) (e* d3 d3))))

(deftest mm-solve-lu (linear-algebra-suite)
  (let* ((a (mx t
              (1 2)
              (3 4)))
         (x (vec 'lla-double 5 6))
         (x-matrix (aops:reshape-col x))
         (b (mm a x))
         (b-matrix (mm a x-matrix))
         (a-lu (lu a))
         (x-solve-lu ))
    (assert-equality #'num= b (vec 'lla-double 17d0 39d0))
    (assert-equality #'num= (solve a b) x)
    (assert-equality #'num= (solve a-lu b) x)
    ;; (assert-equality #'num= (solve a b-matrix) x-matrix)
    ;; (assert-equality #'num= (solve a-lu b-matrix) x-matrix)
    ))

(deftest mm-hermitian (linear-algebra-suite)
  (let* ((a (mx 'lla-double
              (1 2)
              (3 4)))
         (aat (mm a t))
         (ata (mm t a)))
    (assert-equality #'num= aat (hermitian-matrix (mm a (transpose a))))
    (assert-equality #'num= ata (hermitian-matrix (mm (transpose a) a)))))

;; (deftest (linear-algebra-suite)
;;   vector-mm-outer
;;   (let ((a (clo 'lla-double 2 3 5))
;;         (b (clo 'lla-double 7 11))
;;         (c (clo :complex-single #C(1 2) #C(3 5)))
;;         (cc (clo :complex-single :hermitian
;;                  5 #C(13 1) :/
;;                  * 34)))
;;     ;; ;; outer
;;     ;; (assert-equality #'num= (outer a t 0.5)
;;     ;;              (clo 'lla-double :hermitian
;;     ;;                   2 3 5 :/
;;     ;;                   * 4.5 7.5
;;     ;;                   * * 12.5))
;;     ;; (assert-equality #'num= (outer a b 3)
;;     ;;              (clo 'lla-double
;;     ;;                   42 66 :/
;;     ;;                   63 99
;;     ;;                   105 165))
;;     ;; (assert-equality #'num= (outer c t) cc)
;;     ;; (assert-equality #'num= (outer t c) (conjugate-transpose cc))
;;     ;; ;; dot
;;     ;; (assert-equality #'num= (dot b b) 170d0 :test #'=)
;;     ;; (assert-equality #'num= (dot c c) 39d0 :test #'=)
;; ))


(deftest mm-solve-diagonal (linear-algebra-suite)
  (let* ((a (mx 'double-float
              (1 2 3)
              (4 5 6)))
         (b (vec 'double-float 9d0 11))
         (left-d (diagonal-mx 'double-float 7 8))
         (right-d (diagonal-mx 'double-float 9 10 11))
         (left-a (mm left-d a))
         (a-right (mm a right-d)))
    (assert-equality #'num= (mm (aops:as-array left-d) a) left-a)
    (assert-equality #'num= (mm a (aops:as-array right-d)) a-right)
    (assert-equality #'num= (solve left-d left-a) a)
    (assert-equality #'num= (solve left-d (mm left-d b)) b)))

(deftest mm-solve-triangular (linear-algebra-suite)
  (let+ ((u (upper-triangular-mx t
              (1 2)
              (0 3)))
         (l (lower-triangular-mx t
              (1 0)
              (2 3)))
         (b (mx 'lla-double
              (5 6)
              (7 8)))
         (x-u (solve u b))
         (x-l (solve l b)))
    (assert-equality #'num= (mm u x-u) b)
    (assert-equality #'num= (mm l x-l) b)))

;; ;; (defmacro test-update-hermitian% (a x alpha)
;; ;;   (once-only (a x alpha)
;; ;;     `(assert-equality #'num= (update-hermitian ,a ,x ,alpha)
;; ;;                   (e+ ,a (outer ,x t ,alpha))
;; ;;                   :test #'num=)))

;; ;; (deftest (linear-algebra-suite)
;; ;;   update-hermitian
;; ;;   (let ((ad (clo 'lla-double :hermitian
;; ;;                  1 2 :/
;; ;;                  * 4))
;; ;;         (bcd (clo 'lla-complex-double :hermitian
;; ;;                   1 #C(3 4) :/
;; ;;                   * -5))
;; ;;         (xs (clo :single 7 9)))
;; ;;     (test-update-hermitian% ad xs 17)
;; ;;     (test-update-hermitian% bcd xs 17)
;; ;;     (assert-condition error (update-hermitian ad #(1 2 3) 17))
;; ;;     (assert-condition error (update-hermitian ad xs #C(17 9)))))

;; ;; (defmacro test-update-hermitian2% (a x y alpha)
;; ;;   (once-only (a x y alpha)
;; ;;     (with-unique-names (x* y* r1 r2)
;; ;;       `(let* ((,x* (as-column ,x))
;; ;;               (,y* (as-column ,y))
;; ;;               (,r1 (update-hermitian2 ,a ,x ,y ,alpha))
;; ;;               (,r2 (copy-matrix
;; ;;                        (e+ ,a
;; ;;                            (mm ,x* (conjugate-transpose ,y*) ,alpha)
;; ;;                            (mm ,y* (conjugate-transpose ,x*) (conjugate ,alpha)))
;; ;;                        :kind :hermitian)))
;; ;;          (assert-equality #'num= ,r1 ,r2 :test #'num=)))))

;; ;; (deftest (linear-algebra-suite)
;; ;;   update-hermitian2
;; ;;   (let ((ad (clo 'lla-double :hermitian
;; ;;                  1 2 :/
;; ;;                  * 4))
;; ;;         (bcd (clo 'lla-complex-double :hermitian
;; ;;                   1 #C(3 4) :/
;; ;;                   * -5))
;; ;;         (xs (clo :single 7 9))
;; ;;         (yd (clo 'lla-double 5 19)))
;; ;;     (test-update-hermitian2% ad xs yd 17)
;; ;;     (test-update-hermitian2% bcd xs yd 17)
;; ;;     (assert-condition error (update-hermitian2 ad xs #(1 2 3) 17))))

(deftest invert (linear-algebra-suite)
  (let ((m (mx 'lla-double
             (1 2)
             (3 4))))
    (assert-equality #'num= (invert m) ; also tests (invert lu)
        (mx 'lla-double
          (-2 1)
          (1.5 -0.5)))
    (assert-equality #'num= (invert (upper-triangular-matrix m))
        (upper-triangular-mx 'lla-double
          (1 -0.5)
          (0 0.25)))
    (assert-equality #'num= (invert (lower-triangular-matrix m))
        (lower-triangular-mx 'lla-double
          (1 0)
          (-0.75 0.25)))))

;; ;; (deftest (linear-algebra-suite)
;; ;;   eigen
;; ;;   (bind ((a (clo 'lla-double
;; ;;                  1 2 :/
;; ;;                  3 4))
;; ;;          ((:values eigenvalues eigenvectors)
;; ;;           (eigen a :vectors? t :check-real? t))
;; ;;          ((:values ev order) (sort-order eigenvalues #'<)))
;; ;;     ;; we order by eigenvector magnitude
;; ;;     (assert-equality #'num= eigenvalues
;; ;;                  (clo 'lla-double -0.3722813 5.3722813))
;; ;;     (assert-equality #'num= (sub eigenvectors t order)
;; ;;                  (clo 'lla-double
;; ;;                       -0.8245648 -0.4159736 :/
;; ;;                       0.5657675 -0.9093767))))

;; ;; (deftest (linear-algebra-suite)
;; ;;   eigen-single
;; ;;   ;; necessary to test for single type because eigen is so hairy
;; ;;   (bind ((a (clo :single
;; ;;                  1 2 :/
;; ;;                  3 4))
;; ;;          ((:values eigenvalues eigenvectors)
;; ;;           (eigen a :vectors? t :check-real? nil))
;; ;;          ((:values ev order) (sort-order eigenvalues #'<
;; ;;                                          :key #'abs)))
;; ;;     ;; we order by eigenvector magnitude
;; ;;     (assert-equality #'num= eigenvalues
;; ;;                  (clo :complex-single -0.3722813 5.3722813))
;; ;;     (assert-equality #'num= (sub eigenvectors t order)
;; ;;                  (clo :complex-single
;; ;;                       -0.8245648 -0.4159736 :/
;; ;;                       0.5657675 -0.9093767))))

;; ;; (deftest (linear-algebra-suite)
;; ;;   eigen-complex
;; ;;   (bind ((a (clo 'lla-complex-double
;; ;;                  #C(1 2) #C(3 4) :/
;; ;;                  #C(5 6) #C(7 8)))
;; ;;          ((:values eigenvalues eigenvectors)
;; ;;           (eigen a :vectors? t :check-real? nil))
;; ;;          ((:values ev order) (sort-order eigenvalues #'<
;; ;;                                          :key #'abs)))
;; ;;     ;; we order by eigenvector magnitude
;; ;;     (assert-equality #'num= eigenvalues
;; ;;                  (clo 'lla-complex-double
;; ;;                       #C(-0.8845984 -0.7323034)
;; ;;                       #C(8.8845984 10.7323034)))
;; ;;     (assert-equality #'num= (sub eigenvectors t order)
;; ;;                  (clo 'lla-complex-double
;; ;;                       #C(0.8331344 0.0000000) #C(0.3895109 0.0355148) :/
;; ;;                       #C( -0.5526350 -0.0219453) #C(0.9203369 0.0000000)))))

;; ;; (deftest (linear-algebra-suite)
;; ;;   hermitian
;; ;;   (bind ((a (clo 'lla-double
;; ;;                  1 2 :/
;; ;;                  3 4))
;; ;;          (aa (mm t a))
;; ;;          ((:values eigenvalues eigenvectors)
;; ;;           (eigen aa :vectors? t))
;; ;;          ((:values ev order) (sort-order eigenvalues #'<)))
;; ;;     ;; we order by eigenvector magnitude
;; ;;     (assert-equality #'num= ev
;; ;;                  (clo 'lla-double 0.1339313 29.8660687))
;; ;;     (assert-equality #'num= (sub eigenvectors t order)
;; ;;                  (clo 'lla-double
;; ;;                       -0.8174156 0.5760484 :/
;; ;;                       0.5760484 0.8174156))))

(deftest least-squares (linear-algebra-suite)
  (let+ ((x (mx 'lla-double
              (23 23)
              (22 21)
              (25 20)
              (29 32)
              (24 29)))
         (y (vec 'lla-double 67 63 65 94 84))
         (beta (vec 'lla-double 0.7633278 2.2350028))
         (ss 6.724986d0)
         ;; ((:values beta1 ss1 nu1 other1) (least-squares y x :method :svd-d))
         ((&values beta2 ss2 nu2 qr) (least-squares y x :method :qr))
         (raw-var (invert-xx qr))
         (variance (e* (aops:as-array raw-var) (/ ss2 nu2)))
         (*num=-tolerance* 1e-5))
    ;; (assert-equality #'num= beta1 beta)
    (assert-equality #'num= beta2 beta)
    ;; (assert-equality #'num= ss1 ss :test #'num=)
    (assert-equality #'num= ss2 ss)
    ;; (assert-equality #'num= nu1 3 :test #'=)
    (assert-equality #'num= nu2 3 :test #'=)
    ;; (assert-equality #'num= variance (clo 'lla-double :hermitian
    ;;                            0.04035491 * :/
    ;;                            -0.03885797 0.03810950))
    ))

;; ;; (deftest (linear-algebra-suite)
;; ;;   constrained-least-squares
;; ;;   ;;  Taken from the LAPACK documentation
;; ;;   (bind ((x (clo 1 1 1 1 :/
;; ;;                  1 3 1 1
;; ;;                  1 -1 3 1
;; ;;                  1 1 1 3
;; ;;                  1 1 1 -1))
;; ;;          (y (clo 2 1 6 3 1))
;; ;;          (z (clo 1 1 1 -1 :/
;; ;;                  1 -1 1 1
;; ;;                  1 1 -1 1))
;; ;;          (w (clo 1 3 -1))
;; ;;          (*num=-tolerance* 1e-5))
;; ;;     (assert-equality #'num= (constrained-least-squares y x z w)
;; ;;                  (clo :single 0.5 -0.5 1.5 0.5))))

(deftest cholesky (linear-algebra-suite)
  ;; also tests hermitian factorizations
  (let* ((a (hermitian-mx 'lla-double
              (2 -1 0)
              (-1 2 -1)
              (0 -1 2)))
         (l (lower-triangular-mx 'lla-double
              (1.414214)
              (-0.7071068 1.2247449)
              (0.000000 -0.8164966 1.1547005)))
         (c (cholesky a))
         (b (vec 'lla-double 5 7 13))
         (a\b (solve (aops:as-array a) b))
         (a\1 (hermitian-matrix (invert (aops:as-array a))))
         (*lift-equality-test* #'num=))
    (assert-equality #'num= (left-square-root c) l)
    (assert-equality #'num= (solve a b) a\b)
    (assert-equality #'num= (solve c b) a\b)
    (assert-equality #'num= (invert a) a\1)
    (assert-equality #'num= (invert c) a\1)))

(deftest spectral-factorization (linear-algebra-suite)
  (let+ ((a (mm t (mx 'lla-double
                    (1 2)
                    (3 4))))
         (w-true (diagonal-mx 'lla-double 0.1339313 29.8660687))
         (z-true (mx 'lla-double
                   (-0.8174156 0.5760484)
                   (0.5760484  0.8174156)))
         ((&structure-r/o spectral-factorization- z w)
          (spectral-factorization a)))
    (assert-equality #'num= w w-true)
    (assert-equality #'num= z z-true)))

(deftest svd (linear-algebra-suite)
  (let+ ((a (mx 'lla-double
              (0 1)
              (2 3)
              (4 5)))
         (d (diagonal-mx 'lla-double 7.38648 0.66324))
         (u (mx 'lla-double
              (-0.10819  0.90644)
              (-0.48734  0.30958)
              (-0.86649 -0.28729)))
         (v (mx 'lla-double
              (-0.60118 -0.79911)
              (-0.79911  0.60118)))
         (svd1 (svd a))
         (svd2 (svd a :all))
         ((&flet svd-rec (a &optional (vectors :thin))
            (as-array (svd a vectors)))))
    (assert-equality #'num= (svd-d svd1) d)
    (assert-equality #'num= (svd-d svd2) d)
    (assert-equality #'num= (slice (svd-u svd2) t (cons 0 2)) u)
    (assert-equality #'num= (svd-vt svd2) (transpose v))
    (loop repeat 100 do
             (let* ((a0 (+ 2 (random 3)))
                    (a1 (+ 2 (random 3)))
                    (a (aops:generate* 'double-float
                                       (lambda () (random 1d0))
                                       (list a0 a1))))
               (assert-equality #'num= (aops:as-array (svd a :thin)) a)
               (assert-equality #'num= (aops:as-array (svd a :all)) a)))))

;; ;; (deftest (linear-algebra-suite)
;; ;;   svd-rectangular
;; ;;   (bind ((a (clo 'lla-double
;; ;;                  1 2 :/
;; ;;                  3 4
;; ;;                  5 6))
;; ;;          ((:values d u vt) (svd a :left :singular :right :singular)))
;; ;;     (assert-equality #'num= d (clo :diagonal 'lla-double 9.5255181 0.5143006))
;; ;;     (assert-equality #'num= u (clo 'lla-double
;; ;;                         0.2298477  0.8834610 :/
;; ;;                         0.5247448  0.2407825
;; ;;                         0.8196419 -0.4018960))
;; ;;     (assert-equality #'num= vt (clo 'lla-double
;; ;;                          0.6196295 0.7848945 :/
;; ;;                          -0.7848945 0.6196295))))

(deftest tr (linear-algebra-suite)
  (assert-equality #'num= (tr (mx 'lla-double
                                (1 2)
                                (3 4))) 5d0)
  (assert-equality #'num= (tr (hermitian-mx 'lla-double
                                (1 2)
                                (0 9))) 10d0)
  (assert-equality #'num= (tr (mx 'lla-double
                                (1 2)
                                (3 4))) 5d0)
  (assert-equality #'num= (tr (upper-triangular-mx 'lla-double
                                (1 2)
                                (3 4))) 5d0)
  (assert-equality #'num= (tr (diagonal-mx 'lla-double 2 15)) 17d0))

;; (deftest (linear-algebra-suite)
;;   rank
;;   (let ((*lift-equality-test* #'=))
;;     (assert-equality #'num= (rank (clo 1 1 :/ 1 1)) 1)
;;     (assert-equality #'num= (rank (clo 2 4 1 3 :/
;;                             -1 -2 1 0
;;                             0 0 2 2
;;                             3 6 2 5)) 2)))

;; (deftest (linear-algebra-suite)
;;   det
;;   (let ((*num=-tolerance* 1e-4)
;;         (*lift-equality-test* #'num=))
;;     ;; dense
;;     (assert-equality #'num= (det (clo 'lla-double
;;                            1 2 :/
;;                            3 4))
;;                  -2)
;;     ;; upper
;;     (assert-equality #'num= (det (clo :upper
;;                            1 2 :/
;;                            0 4))
;;                  4)
;;     (assert-equality #'num= (det (clo :upper
;;                            1 2 :/
;;                            0 -4))
;;                  -4)
;;     (assert-equality #'num= (det (clo :upper
;;                            1 2 :/
;;                            0 0))
;;                  0)
;;     ;; lower-triangular-mx
;;     (assert-equality #'num= (det (clo :lower
;;                            7 0 :/
;;                            8 12))
;;                  (* 7 12))
;;     (assert-equality #'num= (det (clo :lower
;;                            -7 0 :/
;;                            8 12))
;;                  (* -7 12))
;;     ;; hermitian
;;     ;; (assert-equality #'num= (det (mm t (clo 1 2 :/
;;     ;;                              3 4)))
;;     ;;              4))
;;     ))
