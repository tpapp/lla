(in-package :lla)

;; do vector operations
;; elementwise, element & scalar
;; dot product
;; matrix operations just ride on this, w/ checking conformability
;; names like m+, m-, m*, m/, dot-product
;; NO R-like vector recycling, do that with xref
;; if no type: fall back on xref
;; equality
;; sum
;; product
;; or do a general reduce
;; some of this sould really belong to xref, not lla

;;;; some operations

(defun elements-type (lla-type)
  "Element type of the ELEMENTS a numeric-vector."
  `(simple-array ,(lla-type->lisp-type lla-type) (*)))

(defun elementwise-operation% (op a b &optional integer-to-float-p)
  (let* ((result-type (common-target-type a b))
         (a-elements (elements a))
         (b-elements (elements b))
         (convert-p (and (eq result-type :integer) integer-to-float-p))
         (result-type (if convert-p
                          (if *force-double* :double :single)
                          result-type))
         (lisp-type (lla-type->lisp-type result-type))
         (function (if convert-p
                       (lambda (a b) (coerce (funcall op a b) lisp-type))
                       op)))
    (cons
      (map (elements-type result-type) function a-elements b-elements)
      result-type)))

(defun element-scalar-operation% (op a scalar &key integer-to-float-p flip-p)
  (let* ((scalar-type (lisp-type->lla-type (type-of scalar)))
         (result-type (common-target-type a scalar-type))
         (convert-p (and (eq result-type :integer) integer-to-float-p))
         (result-type (if convert-p
                          (if *force-double* :double :single)
                          result-type))
         (a-elements (elements a))
         (scalar (coerce* scalar result-type))
         (function (if flip-p
                       (lambda (a) (funcall op scalar a))
                       (lambda (a) (funcall op a scalar)))))
    (cons
      (map (elements-type result-type) function a-elements)
      result-type)))

(defgeneric operation% (op a b)
  (:documentation "Elementwise operation on two numeric vectors, or a
  numeric vector and a scalar.  Return the resulting ELEMENTS and
  LLA-TYPE as values.")
  ;; standard cases
  (:method (op (a numeric-vector) (b numeric-vector))
    (elementwise-operation% op a b))
  (:method (op (a numeric-vector) (b number))
    (element-scalar-operation% op a b))
  (:method (op (a number) (b numeric-vector))
    (element-scalar-operation% op b a :flip-p t))
  ;; integers are not closed under /
  (:method ((op (eql #'/)) (a numeric-vector) (b numeric-vector))
    (elementwise-operation% op a b t))
  (:method ((op (eql #'/)) (a numeric-vector) (b number))
    (element-scalar-operation% op a b :integer-to-float-p t))
  (:method ((op (eql #'/)) (a number) (b numeric-vector))
    (element-scalar-operation% op b a :integer-to-float-p t :flip-p t)))

(defmacro define-nv-binary-operation% (op &key integer-to-float-p
                                       (method-name (make-symbol* 'X op)))
  (check-type integer-to-float-p boolean)
  `(flet ((->nv (elements-type-pair)
            (make-nv* (cdr elements-type-pair) (car elements-type-pair))))
     (defmethod ,method-name ((a numeric-vector) (b numeric-vector)
                              &key &allow-other-keys)
       (assert (= (length (elements a)) (length (elements b))) ()
               "A and B don't have equal length.")
       (->nv (elementwise-operation% #',op a b ,integer-to-float-p)))
     (defmethod ,method-name ((a numeric-vector) (b number)
                              &key &allow-other-keys)
       (->nv (element-scalar-operation% #',op a b
                                        :integer-to-float-p ,integer-to-float-p)))
     (defmethod ,method-name ((a number) (b numeric-vector) 
                              &key &allow-other-keys)
       (->nv (element-scalar-operation% #',op b a
                                        :integer-to-float-p ,integer-to-float-p
                                        :flip-p t)))))

(define-nv-binary-operation% +)
(define-nv-binary-operation% -)
(define-nv-binary-operation% *)
(define-nv-binary-operation% / :integer-to-float-p t)


(defmacro define-matrix-binary-operation% (op &key integer-to-float-p 
                                           (method-name (make-symbol* 'X op))
                                           preserve-kind-p 
                                           inherit-zeros-p)
  ;; This is quite a dirty function.  An attempt is made to
  ;; preserve/deduce matrix kind for operations that allow it.  In
  ;; particular:
  ;;
  ;; - Matrix-scalar operations:
  ;;   * and / ALWAYS preserve matrix kind (dense, lower/upper, hermitian)
  ;;   + and - NEVER preserve matrix kind, result is dense
  ;;
  ;; - Elementwise matrix operations:
  ;;   if matrices are of the same kind, so is the result
  ;;   if they are different kind, the result is dense, except for
  ;;   multiplication by lower/upper triangular matrices
  ;;
  ;; - Theoretically, a lower triangular matrix multiplied by an upper
  ;;   triangular one would be diagonal, but here constrain ourselves
  ;;   to returning a dense-matrix-like object.  ??? Maybe we could
  ;;   return a diagonal here, I haven't fixed this in xarray semantics.
  ;;
  ;; !! Currently, we just call SET-RESTRICTED indiscriminately, but
  ;; !! it is not always needed.  Not a main concern at the moment,
  ;; !! treating all the cases right would be a nest of nasty bugs.
  "Define a binary operation OP on two arguments, of which one or two
are matrices.  INTEGER-TO-FLOAT-P should be true iff integers are not
closed under OP.  PRESERVE-KIND-P indicates whether the elementwise
matrix-scalar operations preserves the matrix kind (eg multiplying a
hermitian matrix by a scalar results in a scalar matrix).
INHERIT-ZEROS-P should be true when (= (OP X 0) 0) for every possible
X."
  (check-type integer-to-float-p boolean)
  (check-type preserve-kind-p boolean)
  (check-type inherit-zeros-p boolean)
  `(flet ((->similar-matrix (matrix elements-type-pair kind)
            "Create a matrix similar to MATRIX from ELEMENTS-TYPE-PAIR."
            (make-matrix* (cdr elements-type-pair) (nrow matrix) (ncol matrix)
                          (car elements-type-pair) :kind kind)))
     (defmethod ,method-name ((a dense-matrix-like) (b dense-matrix-like)
                              &key &allow-other-keys)
       (assert (equalp (xdims a) (xdims b)) () "A and B don't have equal dimensions.")
       (set-restricted a)
       (set-restricted b)
       (let* ((a-kind (matrix-kind a))
              (b-kind (matrix-kind b))
              (result-kind 
               (if (eq a-kind b-kind)
                   a-kind
                   ,(if inherit-zeros-p
                        `(cond 
                           ;; if one is lower, other is upper, result
                           ;; will inherit a's kind
                           ((member a-kind '(:lower :upper)) a-kind)
                           ((member b-kind '(:lower :upper)) b-kind)
                           ;; one is Hermitian, other dense
                           (t :dense))
                        :dense))))
         (->similar-matrix a (elementwise-operation% #',op a b
                                                     ,integer-to-float-p)
                           result-kind)))
     (defmethod ,method-name ((a dense-matrix-like) (b number)
                              &key &allow-other-keys)
       (set-restricted a)
       (->similar-matrix a (element-scalar-operation% 
                            #',op a b
                            :integer-to-float-p ,integer-to-float-p)
                         ,(if preserve-kind-p
                              '(matrix-kind a)
                              :dense)))
     (defmethod ,method-name ((a number) (b dense-matrix-like)
                              &key &allow-other-keys)
       (set-restricted b)
       (->similar-matrix b (element-scalar-operation% 
                            #',op b a
                            :integer-to-float-p ,integer-to-float-p
                            :flip-p t)
                         ,(if preserve-kind-p
                              '(matrix-kind b)
                              :dense)))))

(define-matrix-binary-operation% +)
(define-matrix-binary-operation% -)
(define-matrix-binary-operation% * :preserve-kind-p t 
  :inherit-zeros-p t)
(define-matrix-binary-operation% / :preserve-kind-p t)
