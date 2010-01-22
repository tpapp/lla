(in-package :lla)

;;;; transpose
;;;
;;; Methods should take care of returning the correct result type (eg
;;; lower-triangular into upper-triangular, etc).  Helper function
;;; transpose% can be used for implementation.

(defgeneric transpose% (matrix transposed-matrix-kind conjugate-p)
  (:documentation "Return the transpose of MATRIX, which will be of
kind TRANSPOSED-MATRIX-KIND.  SET-RESTRICTED is *NOT* enforced (so the
caller has to decide whether to enforce it).  Meant to be used as a
helper function, *NOT EXPORTED*."))
(expand-for-lla-types (lla-type)
  (let ((complex-p (lla-complex-p lla-type)))
    `(defmethod transpose% ((matrix ,(nv-class lla-type)) transposed-matrix-kind conjugate-p)
       ,@(unless complex-p
           '((declare (ignore conjugate-p))))
       (declare (optimize (speed 3) (safety 0)))
       (bind (((:slots-read-only nrow ncol elements) matrix)
              (transposed (make-nv-elements (length elements) ,lla-type)))
         (declare (fixnum nrow ncol)
                  (type ,(nv-array-type lla-type) elements transposed))
         (dotimes (col ncol)
           (declare (fixnum col))
           (iter
             (declare (iterate:declare-variables))
             (for (the fixnum elements-i) :from (cm-index2 nrow 0 col)
                  :below (cm-index2 nrow nrow col))
             (for (the fixnum transposed-i) :from (cm-index2 ncol col 0) :by ncol)
             (setf (aref transposed transposed-i)
                   ,(if complex-p
                        '(let ((elt (aref elements elements-i)))
                          (if conjugate-p (conjugate elt) elt))
                        '(aref elements elements-i)))))
         (make-matrix* ,lla-type ncol nrow transposed :kind transposed-matrix-kind)))))

(defgeneric transpose (matrix &optional conjugate-p)
  (:documentation "Return the transpose of a matrix.")
  (:method ((matrix dense-matrix) &optional (conjugate-p t))
    (transpose% matrix :dense conjugate-p))
  (:method ((matrix upper-triangular-matrix) &optional (conjugate-p t))
    (transpose% matrix :lower-triangular conjugate-p))
  (:method ((matrix lower-triangular-matrix) &optional (conjugate-p t))
    (transpose% matrix :upper-triangular conjugate-p))
  (:method ((matrix hermitian-matrix) &optional (conjugate-p t))
    (set-restricted matrix)
    (case (lla-type matrix)
      ((:complex-single :complex-double)
         (transpose% matrix :hermitian conjugate-p))
      (otherwise
         (copy-matrix matrix :copy-p t)))))


;; (defmethod x* ((a dense-matrix-like) (b number) &key)
;;   ;;; !!!! this assumes that all restricted elements are zeros
;;   (if (typep a 'restricted-elements)
;;       (make-instance (class-of a) :nrow (nrow a) :ncol (ncol a)
;;                      :data (x* (data a) b) :restricted-set-p (restricted-set-p a))
;;       (make-instance (class-of a) :nrow (nrow a) :ncol (ncol a)
;;                      :data (x* (data a) b))))

;; (defmethod x* ((a number) (b dense-matrix-like) &key)
;;   (x* b a))


;;;; extraction methods for factorization components

(defgeneric matrix-from-first-rows (vector lla-type m nrhs n &optional kind)
  (:documentation "Extract & return (as KIND, default is :DENSE) the
first N rows of an MxNRHS matrix, given as a Lisp vector with
column-major indexing.  NOTE: needed to interface to LAPACK routines
like xGELS."))
(expand-for-lla-types (lla-type)
  (let ((array-type (nv-array-type lla-type)))
    `(defmethod matrix-from-first-rows (vector (lla-type (eql ,lla-type)) m nrhs n
                                        &optional (kind :dense))
       ;; It is assumed that NV's ELEMENTS has the correct type.
       (declare (optimize speed (safety 0))
                (type ,array-type vector)
                (type fixnum m nrhs n))
       (let ((result (make-nv-elements (the fixnum (* n nrhs)) ,lla-type)))
         (dotimes (col nrhs)
           (declare (fixnum col))
           (iter
             (declare (iterate:declare-variables))
             (repeat n)
             (for (the fixnum vector-index) :from (the fixnum (* col m)))
             (for (the fixnum result-index) :from (the fixnum (* col n)))
             (setf (aref result result-index) (aref vector vector-index))))
         (make-matrix* lla-type n nrhs result :kind kind)))))

(defmethod factorization-component ((mf qr) (component (eql :R)) &key copy-p)
  (declare (ignore copy-p))
  (bind (((:slots-read-only qr-matrix) mf)
         ((:slots-read-only nrow ncol elements) qr-matrix))
    (matrix-from-first-rows elements (lla-type qr-matrix) nrow ncol ncol :upper-triangular)))

(defmethod factorization-component ((mf cholesky) component &key (copy-p nil))
  (flet ((copy-maybe (matrix)
           (if copy-p
               (copy-matrix matrix :copy-p t)
               matrix)))
    (bind (((:slots-read-only factor) mf))
      (etypecase factor
        (lower-triangular-matrix
           (ecase component
             (:U (copy-maybe (transpose factor)))
             (:L factor)))
        (upper-triangular-matrix
           (ecase component
             (:U factor)
             (:L (copy-maybe (transpose factor)))))))))

;;;; elementwise and scalar operations

(defmacro define-matrix-operations (op &optional (method-name (make-symbol* 'X op)))
  "Define bivariate (both elementwise and with scalar) operations for
matrices.  Note: just dispatches to vector operations."
  `(progn
     (defmethod ,method-name ((a dense-matrix-like) (b dense-matrix-like)
                              &key element-type)
       (declare (ignore element-type))
       (set-restricted a)
       (set-restricted b)
       (bind (((:slots-read-only nrow ncol) a))
         (assert (and (= nrow (nrow b)) (= ncol (ncol b))) ()
                 "Incompatible dimensions.")
         (vector->matrix (,method-name (copy-nv a) (copy-nv b)) nrow ncol)))
     (defmethod ,method-name ((a dense-matrix-like) (b number)
                              &key element-type)
       (declare (ignore element-type))
       (set-restricted a)
       (bind (((:slots-read-only nrow ncol) a))
         (vector->matrix (,method-name (copy-nv a) b) nrow ncol)))
     (defmethod ,method-name ((a number) (b dense-matrix-like)
                              &key element-type)
       (declare (ignore element-type))
       (set-restricted b)
       (bind (((:slots-read-only nrow ncol) b))
         (vector->matrix (,method-name a (copy-nv b)) nrow ncol)))))

(define-matrix-operations +)
(define-matrix-operations *)
(define-matrix-operations -)
(define-matrix-operations /)

     
;;;; stacking

(defun stack-vertically (&rest arguments)
  "Stack arguments vertically, converting to a common type.  A vector
  is interpreted as a row matrix."
  (let* ((matrices (mapcar (lambda (arg)
                             (etypecase arg
                               (diagonal (diagonal->matrix arg))
                               (dense-matrix-like (set-restricted arg))
                               (numeric-vector (vector->row arg))))
                           arguments))
         (common-type (apply #'common-target-type (mapcar #'lla-type matrices)))
         (nrow-total (reduce #'+ matrices :key #'nrow))
         (ncol (ncol (first matrices)))
         (result (make-matrix nrow-total ncol common-type))
         (result-elements (elements result))
         (row 0))
    (iter
      (for matrix :in matrices)
      (for nrow := (nrow matrix))
      (for elements := (elements matrix))
      (for lla-type := (lla-type matrix))
      (assert (= ncol (ncol matrix)) () "Matrix columns (or vector lengths) do not match.")
      (iter
        (for col :from 0 :below ncol)
        (copy-elements-into elements lla-type (cm-index2 nrow 0 col)
                            result-elements common-type (cm-index2 nrow-total row col)
                            nrow))
      (incf row nrow))
    result))

(defun stack-horizontally (&rest arguments)
  (declare (ignore arguments))
  (error "This function needs to be implemented."))
