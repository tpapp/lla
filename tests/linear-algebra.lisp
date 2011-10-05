;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla-tests)

(deftestsuite linear-algebra-tests (lla-tests)
  ()
  (:equality-test #'==))

;;;; linear algebra

(addtest (linear-algebra-tests)
  mm
  (let ((a (dense 'lla-double
             (1 2)
             (3 4)
             (5 6)))
        (b2 (vec 'lla-double 1 2))
        (b3 (vec 'lla-double 1 2 3))
        (i2 (dense 'lla-double 
              (1 0)
              (0 1)))
        (i3 (dense 'lla-double
              (1 0 0)
              (0 1 0)
              (0 0 1))))
    (ensure-same (mm a i2) a)
    (ensure-same (mm i3 a) a)
    (ensure-same (mm a b2) (vec 'lla-double 5 11 17))
    (ensure-same (mm b3 a) (vec 'lla-double 22 28))
    (ensure-error (mm a b3))
    ;; (ensure-error (mm b3 b3))
    (ensure-error (mm b2 b3))))

(addtest (linear-algebra-tests)
  mm-dot
  (let+ ((a (vec 'lla-double 2 3 5))
         (b (vec 'lla-complex-double 1 #C(2 1) 3))
         ((&flet mm-vec (a b)
            (sum (map 'vector (lambda (a b)
                                (* a (conjugate b)))
                      a b))))
         ((&macrolet test (a1 b1 a2 b2)
            `(ensure-same (mm ,a1 ,b1)
                          (mm-vec ,a2 ,b2)))))
    (test a b a b)
    (test b a b a)
    (test a t a a)
    (test b t b b)
    (test t a a a)
    (test t b b b)))

(addtest (linear-algebra-tests)
  outer
  (let+ ((a (vec 'lla-double 2 3))
         (b (vec 'lla-complex-double 1 #C(2 1) 9))
         ((&flet outer2 (a b &optional hermitian?)
            (let ((result (mm (as-column a)
                              (as-row (map1 #'conjugate b)))))
              (if hermitian?
                  (convert-matrix 'hermitian result)
                  result))))
         (*lift-equality-test* #'==))
    (ensure-same (outer a b) (outer2 a b))
    (ensure-same (outer a t) (outer2 a a t))
    (ensure-same (outer t a) (outer2 a a t))
    (ensure-same (outer (vec 'lla-complex-double 1 #C(2 1) 9) t)
                 (hermitian 'lla-complex-double
                   (1)
                   (#C(2 1) 5)
                   (9 #C (18 -9) 81)))))

(addtest (linear-algebra-tests)
  mm-diagonal
  (let ((a (dense 'lla-double
             (1 2)
             (3 4)
             (5 6)))
        (d2 (diag 'lla-double 2 3))
        (d3 (diag 'lla-double 2 3 5))
        (d3*a (dense 'lla-double
                (2 4)
                (9 12)
                (25 30)))
        (a*d2 (dense 'lla-double
                (2 6)
                (6 12)
                (10 18)))
        (*lift-equality-test* #'equalp))
    (ensure-same (mm a d2) a*d2)
    (ensure-same (mm d3 a) d3*a)
    (ensure-error (mm a d3))
    (ensure-error (mm d2 a))
    (ensure-same (mm d3 d3) (e* d3 d3))
    (ensure-same (mm d3 t) (e* d3 d3))  
    (ensure-same (mm t d3) (e* d3 d3))))

(addtest (linear-algebra-tests)
  mm-solve-lu
  (let* ((a (dense t
              (1 2)
              (3 4)))
         (x (vec 'lla-double 5 6))
         (x-matrix (as-column x))
         (b (mm a x))
         (b-matrix (mm a x-matrix))
         (a-lu (lu a))
         (x-solve-lu ))
    (ensure-same b (vec 'lla-double 17d0 39d0))
    (ensure-same (solve a b) x)
    (ensure-same (solve a-lu b) x)
    (ensure-same (solve a b-matrix) x-matrix)
    (ensure-same (solve a-lu b-matrix) x-matrix)))

(addtest (linear-algebra-tests)
  mm-hermitian
  (let* ((a (dense 'lla-double
              (1 2)
              (3 4)))
         (aat (mm a t))
         (ata (mm t a))
         (*lift-equality-test* #'==))
    (ensure-same aat (convert-matrix 'hermitian (mm a (transpose a))))
    (ensure-same ata (convert-matrix 'hermitian (mm (transpose a) a)))))

;; (addtest (linear-algebra-tests)
;;   vector-mm-outer
;;   (let ((a (clo 'lla-double 2 3 5))
;;         (b (clo 'lla-double 7 11))
;;         (c (clo :complex-single #C(1 2) #C(3 5)))
;;         (cc (clo :complex-single :hermitian
;;                  5 #C(13 1) :/
;;                  * 34)))
;;     ;; ;; outer
;;     ;; (ensure-same (outer a t 0.5)
;;     ;;              (clo 'lla-double :hermitian
;;     ;;                   2 3 5 :/
;;     ;;                   * 4.5 7.5
;;     ;;                   * * 12.5))
;;     ;; (ensure-same (outer a b 3)
;;     ;;              (clo 'lla-double
;;     ;;                   42 66 :/
;;     ;;                   63 99
;;     ;;                   105 165))
;;     ;; (ensure-same (outer c t) cc)
;;     ;; (ensure-same (outer t c) (conjugate-transpose cc))
;;     ;; ;; dot
;;     ;; (ensure-same (dot b b) 170d0 :test #'=)
;;     ;; (ensure-same (dot c c) 39d0 :test #'=)
;; ))


;; ;; (addtest (linear-algebra-tests)
;; ;;   mm-solve-diagonal
;; ;;   (let* ((a (clo 1 2 3 :/
;; ;;                  4 5 6))
;; ;;          (b (clo 'lla-double 9 11))
;; ;;          (left-d (clo :diagonal 7 8))
;; ;;          (right-d (clo :diagonal 9 10 11))
;; ;;          (left-a (mm left-d a))
;; ;;          (a-right (mm a right-d)))
;; ;;     (ensure-same (mm (as-matrix left-d) a) left-a)
;; ;;     (ensure-same (mm a (as-matrix right-d)) a-right)
;; ;;     (ensure-same (solve left-d left-a) a)
;; ;;     (ensure-same (solve left-d (mm left-d b)) b)))


(addtest (linear-algebra-tests)
  mm-solve-triangular
  (let+ ((u (upper t
              (1 2)
              (0 3)))
         (l (lower t
              (1 0)
              (2 3)))
         (b (dense 'lla-double
              (5 6)
              (7 8)))
         (x-u (solve u b))
         (x-l (solve l b)))
    (ensure-same (mm u x-u) b)
    (ensure-same (mm l x-l) b)))

;; ;; (defmacro test-update-hermitian% (a x alpha)
;; ;;   (once-only (a x alpha)
;; ;;     `(ensure-same (update-hermitian ,a ,x ,alpha)
;; ;;                   (e+ ,a (outer ,x t ,alpha))
;; ;;                   :test #'==)))

;; ;; (addtest (linear-algebra-tests)
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
;; ;;     (ensure-error (update-hermitian ad #(1 2 3) 17))
;; ;;     (ensure-error (update-hermitian ad xs #C(17 9)))))

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
;; ;;          (ensure-same ,r1 ,r2 :test #'==)))))

;; ;; (addtest (linear-algebra-tests)
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
;; ;;     (ensure-error (update-hermitian2 ad xs #(1 2 3) 17))))

(addtest (linear-algebra-tests)
  invert
  (let ((m (dense 'lla-double
             (1 2)
             (3 4))))
    (flet ((invert (kind)
             (invert (convert-matrix kind m :copy? t))))
      (ensure-same (invert 'dense)           ; also tests (invert lu)
                   (dense 'lla-double
                     (-2 1)
                     (1.5 -0.5)))
      (ensure-same (invert 'upper)
                   (upper 'lla-double
                     (1 -0.5)
                     (0 0.25)))
      (ensure-same (invert 'lower)
                   (lower 'lla-double
                     (1 0)
                     (-0.75 0.25))))))

;; ;; (addtest (linear-algebra-tests)
;; ;;   eigen
;; ;;   (bind ((a (clo 'lla-double
;; ;;                  1 2 :/
;; ;;                  3 4))
;; ;;          ((:values eigenvalues eigenvectors)
;; ;;           (eigen a :vectors? t :check-real? t))
;; ;;          ((:values ev order) (sort-order eigenvalues #'<)))
;; ;;     ;; we order by eigenvector magnitude
;; ;;     (ensure-same eigenvalues
;; ;;                  (clo 'lla-double -0.3722813 5.3722813))
;; ;;     (ensure-same (sub eigenvectors t order)
;; ;;                  (clo 'lla-double
;; ;;                       -0.8245648 -0.4159736 :/
;; ;;                       0.5657675 -0.9093767))))

;; ;; (addtest (linear-algebra-tests)
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
;; ;;     (ensure-same eigenvalues
;; ;;                  (clo :complex-single -0.3722813 5.3722813))
;; ;;     (ensure-same (sub eigenvectors t order)
;; ;;                  (clo :complex-single
;; ;;                       -0.8245648 -0.4159736 :/
;; ;;                       0.5657675 -0.9093767))))

;; ;; (addtest (linear-algebra-tests)
;; ;;   eigen-complex
;; ;;   (bind ((a (clo 'lla-complex-double
;; ;;                  #C(1 2) #C(3 4) :/
;; ;;                  #C(5 6) #C(7 8)))
;; ;;          ((:values eigenvalues eigenvectors)
;; ;;           (eigen a :vectors? t :check-real? nil))
;; ;;          ((:values ev order) (sort-order eigenvalues #'<
;; ;;                                          :key #'abs)))
;; ;;     ;; we order by eigenvector magnitude
;; ;;     (ensure-same eigenvalues
;; ;;                  (clo 'lla-complex-double
;; ;;                       #C(-0.8845984 -0.7323034)
;; ;;                       #C(8.8845984 10.7323034)))
;; ;;     (ensure-same (sub eigenvectors t order)
;; ;;                  (clo 'lla-complex-double
;; ;;                       #C(0.8331344 0.0000000) #C(0.3895109 0.0355148) :/
;; ;;                       #C( -0.5526350 -0.0219453) #C(0.9203369 0.0000000)))))

;; ;; (addtest (linear-algebra-tests)
;; ;;   hermitian
;; ;;   (bind ((a (clo 'lla-double
;; ;;                  1 2 :/
;; ;;                  3 4))
;; ;;          (aa (mm t a))
;; ;;          ((:values eigenvalues eigenvectors)
;; ;;           (eigen aa :vectors? t))
;; ;;          ((:values ev order) (sort-order eigenvalues #'<)))
;; ;;     ;; we order by eigenvector magnitude
;; ;;     (ensure-same ev
;; ;;                  (clo 'lla-double 0.1339313 29.8660687))
;; ;;     (ensure-same (sub eigenvectors t order)
;; ;;                  (clo 'lla-double
;; ;;                       -0.8174156 0.5760484 :/
;; ;;                       0.5760484 0.8174156))))

(addtest (linear-algebra-tests)
  least-squares
  (let+ ((x (dense 'lla-double 
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
         (variance (e* (as-array raw-var) (/ ss2 nu2)))
         (*==-tolerance* 1e-5))
    ;; (ensure-same beta1 beta)
    (ensure-same beta2 beta)
    ;; (ensure-same ss1 ss :test #'==)
    (ensure-same ss2 ss)
    ;; (ensure-same nu1 3 :test #'=)
    (ensure-same nu2 3 :test #'=)
    ;; (ensure-same variance (clo 'lla-double :hermitian
    ;;                            0.04035491 * :/
    ;;                            -0.03885797 0.03810950))
    ))

;; ;; (addtest (linear-algebra-tests)
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
;; ;;          (*==-tolerance* 1e-5))
;; ;;     (ensure-same (constrained-least-squares y x z w)
;; ;;                  (clo :single 0.5 -0.5 1.5 0.5))))

(addtest (linear-algebra-tests)
  cholesky                       ; also tests hermitian factorizations
  (let* ((a (hermitian 'lla-double
              (2 -1 0)
              (-1 2 -1)
              (0 -1 2)))
         (l (lower 'lla-double
              (1.414214)
              (-0.7071068 1.2247449)
              (0.000000 -0.8164966 1.1547005)))
         (c (cholesky a))
         (b (vec 'lla-double 5 7 13))
         (a\b (solve (as-array a) b))
         (a\1 (convert-matrix 'hermitian (invert (as-array a))))
         (*lift-equality-test* #'==))
    (ensure-same (left-square-root c) l)
    (ensure-same (solve a b) a\b)
    (ensure-same (solve c b) a\b)
    (ensure-same (invert a) a\1)
    (ensure-same (invert c) a\1)))

(addtest (linear-algebra-tests)
  spectral-factorization
  (let+ ((a (mm t (dense 'lla-double
                    (1 2)
                    (3 4))))
         (w-true (diag 'lla-double 0.1339313 29.8660687))
         (z-true (dense 'lla-double
                   (-0.8174156 0.5760484)
                   (0.5760484  0.8174156)))
         ((&structure-r/o spectral-factorization- z w)
         (spectral-factorization a))
         (*lift-equality-test* #'==))
    (ensure-same w w-true)
    (ensure-same z z-true)))

(addtest (linear-algebra-tests)
  svd
  (let+ ((a (dense 'lla-double
              (0 1)
              (2 3)
              (4 5)))
         (d (diag 'lla-double 7.38648 0.66324))
         (u (dense 'lla-double 
              (-0.10819  0.90644)
              (-0.48734  0.30958)
              (-0.86649 -0.28729)))
         (v (dense 'lla-double
              (-0.60118 -0.79911)
              (-0.79911  0.60118)))
         (svd1 (svd a))
         (svd2 (svd a :all))
         (*lift-equality-test* #'==)
         ((&flet svd-rec (a &optional (vectors :thin))
            (as-array (svd a vectors)))))
    (ensure-same (svd-d svd1) d)
    (ensure-same (svd-d svd2) d)
    (ensure-same (sub (svd-u svd2) t (cons 0 2)) u)
    (ensure-same (svd-vt svd2) (transpose v))
    (loop repeat 100 do
      (let* ((a0 (+ 2 (random 3)))
             (a1 (+ 2 (random 3)))
             (a (filled-array (list a0 a1) (lambda () (random 1d0))
                              'double-float)))
        (ensure-same (as-array (svd a :thin)) a)
        (ensure-same (as-array (svd a :all)) a)))))

;; ;; (addtest (linear-algebra-tests)
;; ;;   svd-rectangular
;; ;;   (bind ((a (clo 'lla-double 
;; ;;                  1 2 :/
;; ;;                  3 4
;; ;;                  5 6))
;; ;;          ((:values d u vt) (svd a :left :singular :right :singular)))
;; ;;     (ensure-same d (clo :diagonal 'lla-double 9.5255181 0.5143006))
;; ;;     (ensure-same u (clo 'lla-double
;; ;;                         0.2298477  0.8834610 :/
;; ;;                         0.5247448  0.2407825
;; ;;                         0.8196419 -0.4018960))
;; ;;     (ensure-same vt (clo 'lla-double
;; ;;                          0.6196295 0.7848945 :/
;; ;;                          -0.7848945 0.6196295))))

(addtest (linear-algebra-tests)
  tr
  (let ((*lift-equality-test* #'=))
   (ensure-same (tr (dense 'lla-double
                      (1 2)
                      (3 4))) 5d0)
    (ensure-same (tr (hermitian 'lla-double
                       (1 2)
                       (0 9))) 10d0)
    (ensure-same (tr (dense 'lla-double
                       (1 2)
                       (3 4))) 5d0)
    (ensure-same (tr (upper 'lla-double
                       (1 2)
                       (3 4))) 5d0)
    (ensure-same (tr (diag 'lla-double 2 15)) 17d0)))

;; (addtest (linear-algebra-tests)
;;   rank
;;   (let ((*lift-equality-test* #'=))
;;     (ensure-same (rank (clo 1 1 :/ 1 1)) 1)
;;     (ensure-same (rank (clo 2 4 1 3 :/
;;                             -1 -2 1 0
;;                             0 0 2 2
;;                             3 6 2 5)) 2)))

;; (addtest (linear-algebra-tests)
;;   det
;;   (let ((*==-tolerance* 1e-4)
;;         (*lift-equality-test* #'==))
;;     ;; dense
;;     (ensure-same (det (clo 'lla-double
;;                            1 2 :/
;;                            3 4))
;;                  -2)
;;     ;; upper
;;     (ensure-same (det (clo :upper
;;                            1 2 :/
;;                            0 4))
;;                  4)
;;     (ensure-same (det (clo :upper
;;                            1 2 :/
;;                            0 -4))
;;                  -4)
;;     (ensure-same (det (clo :upper
;;                            1 2 :/
;;                            0 0))
;;                  0)
;;     ;; lower
;;     (ensure-same (det (clo :lower
;;                            7 0 :/
;;                            8 12))
;;                  (* 7 12))
;;     (ensure-same (det (clo :lower
;;                            -7 0 :/
;;                            8 12))
;;                  (* -7 12))
;;     ;; hermitian
;;     ;; (ensure-same (det (mm t (clo 1 2 :/
;;     ;;                              3 4)))
;;     ;;              4))
;;     ))
