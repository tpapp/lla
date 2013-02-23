(in-package #:lla-tests)

(defsuite blas-suite (tests))

(defun submatrix (matrix height width)
  (let ((a (make-array (list height width)
                       :element-type (array-element-type matrix))))
    (dotimes (row height)
      (dotimes (col width)
        (setf (aref a row col) (aref matrix row col))))
    a))

(defun sloppy-random-array (type height width)
  (let ((a* (apply #'random-array type (mapcar (lambda (dim)
                                                 (+ dim (random 3)))
                                               (list height width)))))
    (values (submatrix a* height width) a*)))

(deftest gemm (blas-suite)
  (loop repeat 100 do
    (let+ ((m (1+ (random 10)))
           (n (1+ (random 10)))
           (k (1+ (random 10)))
           (transpose-a? (if (zerop (random 2)) t nil))
           (transpose-b? (if (zerop (random 2)) t nil))
           ((&values a a*)
            (if transpose-a?
                (sloppy-random-array 'double-float k m)
                (sloppy-random-array 'double-float m k)))
           ((&values b b*)
            (if transpose-b?
                (sloppy-random-array 'double-float n k)
                (sloppy-random-array 'double-float k n)))
           ((&values c c*)
            (sloppy-random-array 'double-float m n))
           (alpha (random 10d0))
           (beta (random 10d0))
           (result (e+ (e* alpha (mm (if transpose-a? (transpose a) a)
                                     (if transpose-b? (transpose b) b)))
                       (e* beta c))))
      (gemm! alpha a* b* beta c*
             :transpose-a? transpose-a?
             :transpose-b? transpose-b?
             :m m :n n :k k)
      (assert-equality #'num= (apply #'submatrix c* (array-dimensions c))
          result))))
