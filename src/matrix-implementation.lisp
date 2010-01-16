(in-package :lla)

;;;; Implementation of basic matrix interface and some methods.
;;;
;;; LAPACK-specific code is in other files.

;;; helper functions

(declaim (inline cm-index2))
(defun cm-index2 (nrow row col)
  "Calculate column-major index, without error checking.  Inlined."
  (the fixnum (+ (the fixnum (* nrow col)) row)))

;;; set-restricted* methods
;;;
;;; !!! write optimized versions

(expand-for-lla-types (lla-type)
  (let ((class (matrix-class :upper-triangular lla-type))
        (zero (coerce 0 (lla-type->lisp-type lla-type))))
    `(defmethod set-restricted* ((matrix ,class))
       (declare (optimize speed))
       (ensure-unshared matrix)
       (bind (((:slots-read-only nrow ncol elements) matrix))
         ;; set the lower triangle (below diagonal) to 0
         (declare (fixnum nrow ncol)
                  (type ,(nv-array-type lla-type) elements))
         (dotimes (col ncol)
           (declare (fixnum col ncol))
           (iter
             (declare (iterate:declare-variables))
             (for (the fixnum index)
                  :from (1+ (cm-index2 nrow col col))
                  :below (cm-index2 nrow nrow col))
             (setf (aref elements index) ,zero))))
       matrix)))

(defmethod set-restricted* ((matrix lower-triangular-matrix))
  (ensure-unshared matrix)
  (bind (((:slots-read-only nrow ncol elements) matrix)
         (zero (coerce 0 (array-element-type elements))))
    ;; set the upper triangle (above diagonal) to 0
    (dotimes (col ncol)
      (iter
        (for index
          :from (cm-index2 nrow 0 col)
          :below (cm-index2 nrow col col))
        (setf (aref elements index) zero))))
  matrix)

(defmethod set-restricted* ((matrix hermitian-matrix))
  (ensure-unshared matrix)
  (bind (((:slots-read-only nrow ncol elements) matrix))
    ;; set the lower triangle (below diagonal) to conjugate of the
    ;; elements in the upper triangle
    (dotimes (col ncol)
      (iter
        (for row :from col :below nrow)
        (for index
          :from (cm-index2 nrow col col)
          :below (cm-index2 nrow nrow col))
        (setf (aref elements index)
              (conjugate (aref elements (cm-index2 nrow col row)))))))
  matrix)


;;;; General XARRAY interface.
;;;
;;; Notes: XELTTYPE is already defined for the NUMERIC-VECTOR
;;; superclasses.

(defmethod xrank ((matrix dense-matrix-like))
  2)

(defmethod xdims ((matrix dense-matrix-like))
  (list (nrow matrix) (ncol matrix)))

(defmethod xdim ((matrix dense-matrix-like) axis-number)
  (ecase axis-number
    (0 (nrow matrix))
    (1 (ncol matrix))))

(defmethod xsize ((matrix dense-matrix-like))
  (* (nrow matrix) (ncol matrix)))

;;;; !!!! slot access could be optimized by dispatching on subclasses
;;;;      according to type
;;;;
;;;;

;;; xref for dense matrices

(declaim (inline matrix-xref% matrix-setf-xref%))
(defun matrix-xref% (matrix subscripts)
  "Reader function for matrices.  Meant to be inlined.  Not exported."
  (bind (((row col) subscripts))
    (with-slots (nrow ncol elements) matrix
      (check-index row nrow)
      (check-index col ncol)
      (aref elements (cm-index2 nrow row col)))))
(defun matrix-setf-xref% (value matrix subscripts)
  "Setter function for matrices.  Meant to be inlined.  Not exported."
  (ensure-unshared matrix)
  (bind (((row col) subscripts))
    (with-slots (nrow ncol elements) matrix
      (check-index row nrow)
      (check-index col ncol)
      (setf (aref elements (cm-index2 nrow row col))
	    value))))

(defmethod xref ((matrix dense-matrix) &rest subscripts)
  (matrix-xref% matrix subscripts))

(defmethod (setf xref) (value (matrix dense-matrix) &rest subscripts)
  (matrix-setf-xref% value matrix subscripts))

;;; xref for upper triangular matrices

(defmethod xref ((matrix upper-triangular-matrix) &rest subscripts)
  (bind (((row col) subscripts))
    (with-slots (nrow ncol elements) matrix
      (check-index row nrow)
      (check-index col ncol)
      (if (<= row col)
          (aref elements (cm-index2 nrow row col))
          (coerce 0 (xelttype matrix))))))

(defmethod (setf xref) (value (matrix upper-triangular-matrix) &rest subscripts)
  (bind (((row col) subscripts))
    (with-slots (nrow ncol elements) matrix
      (check-index row nrow)
      (check-index col ncol)
      (if (<= row col)
          (progn
            (ensure-unshared matrix)
            (setf (aref elements (cm-index2 nrow row col))
                  value))
          (unless (zerop value)
            (error 'xref-setting-readonly))))))

;;; xref for lower triangular matrices

(defmethod xref ((matrix lower-triangular-matrix) &rest subscripts)
  (bind (((row col) subscripts))
    (with-slots (nrow ncol elements) matrix
      (check-index row nrow)
      (check-index col ncol)
      (if (>= row col)
          (aref elements (cm-index2 nrow row col))
          (coerce 0 (xelttype matrix))))))

(defmethod (setf xref) (value (matrix lower-triangular-matrix) &rest subscripts)
  (bind (((row col) subscripts))
    (with-slots (nrow ncol elements) matrix
      (check-index row nrow)
      (check-index col ncol)
      (if (>= row col)
          (progn
            (ensure-unshared matrix)
            (setf (aref elements (cm-index2 nrow row col))
                  value))
          (unless (zerop value)
            (error 'xref-setting-readonly))))))

;;; xref for hermitian matrices

(defmethod xref ((matrix hermitian-matrix) &rest subscripts)
  (bind (((row col) subscripts))
    (with-slots (nrow ncol elements) matrix
      (check-index row nrow)
      (check-index col ncol)
      (if (<= row col)
          (aref elements (cm-index2 nrow row col))
          (conjugate (aref elements (cm-index2 nrow col row)))))))

(defmethod (setf xref) (value (matrix hermitian-matrix) &rest subscripts)
  (ensure-unshared matrix)
  (bind (((row col) subscripts))
    (with-slots (nrow ncol elements restricted-set-p) matrix
      (check-index row nrow)
      (check-index col ncol)
      (setf restricted-set-p nil)
      (if (<= row col)
          (setf (aref elements (cm-index2 nrow row col)) value)
          (setf (aref elements (cm-index2 nrow col row)) (conjugate value))))))

;;;; xref for factorizations

(defmethod xref ((matrix matrix-factorization) &rest subscripts)
  (matrix-xref% matrix subscripts))

(defmethod (setf xref) (value (matrix matrix-factorization) &rest subscripts)
  (matrix-setf-xref% value matrix subscripts))

;;;; matrix creation

(declaim (inline make-matrix*))
(defun make-matrix* (lla-type nrow ncol elements &key (kind :dense)
                     (shared-p nil) (restricted-set-p nil))
  "Create a matrix with given ELEMENTS, TYPE, LLA-TYPE and dimensions.
Note that there is no type checking, and elements are not copied: this
is effectively shorthand for a MAKE-INSTANCE call.  For internal use,
not exported."
  (aprog1 (make-instance (matrix-class kind lla-type) :nrow nrow :ncol ncol
                         :elements elements :shared-p shared-p)
    (set-restricted-set-p it restricted-set-p)))
    

(defun make-matrix (nrow ncol lla-type &key (kind :dense) (initial-element 0))
  "Create a matrix with given parameters, optionally initialized with
INITIAL-ELEMENTs (where it makes sense, for restricted matrices)."
  (make-matrix* lla-type nrow ncol
                (make-nv-elements (* nrow ncol) lla-type initial-element)
                :kind kind))

(defun matrix-elements-from-sequence (ncol lla-type sequence)
  "Return a lisp vector comforming to LLA-TYPE which has initial
contents derived from 2d Lisp ARRAY.  Return (VALUES ELEMENTS NROW).
If NCOL is zero, return a row matrix."
  (bind ((length (length sequence))
         ((:values nrow ncol) (if (zerop ncol)
                                  (values length 1)
                                  (bind (((:values nrow remainder) (floor length ncol)))
                                    (unless (zerop remainder)
                                      (error "Length of sequence (~A) is not ~
                                              a multiple of ncol (~A)." length ncol))
                                    (values nrow ncol))))
	 (lisp-type (lla-type->lisp-type lla-type))
	 (elements (make-nv-elements length lla-type))
         (col 0)
         (row 0))
    (flet ((store-element (x)
             ;; store elements, traversing the elements in a column-major order
             (setf (aref elements (cm-index2 nrow row col)) (coerce x lisp-type))
             (incf col)
             (when (= col ncol)
               (incf row)
               (setf col 0))))
      (etypecase sequence
        (list (dolist (x sequence)
                (store-element x)))
        (vector (iter
                  (for x :in-vector sequence)
                  (store-element x)))))
    (values elements nrow)))

(defun create-matrix (ncol initial-contents &key (kind :dense) lla-type)
  "Create matrix of given TYPE vector with given initial contents (a
sequence).  Unless LLA-TYPE is given, it is inferred from the
elements.  NCOL gives the number of columns, while the number of rows
is inferred from the length of the sequence.  An error is signalled if
there are remainder elements.

Usage note: This is a convenience function for easily creation of
matrices.  Also see *force-float*."
  (bind ((lla-type (infer-lla-type lla-type initial-contents *force-float*))
         ((:values elements nrow)
          (matrix-elements-from-sequence ncol lla-type initial-contents)))
    (make-matrix* lla-type nrow ncol elements :kind kind)))

(defun copy-matrix (matrix &key (kind (matrix-kind matrix))
                    (destination-type (lla-type matrix)) (copy-p nil))
  "Copy or convert matrix to the given kind and destination-type.
Copying is forced when COPY-P."
  (set-restricted matrix)
  (bind (((:slots-read-only nrow ncol) matrix)
         ((:values elements shared-p) (copy-elements* matrix destination-type copy-p)))
    (make-matrix* (lla-type matrix) nrow ncol elements
                  :kind kind :shared-p shared-p
                  :restricted-set-p t)))

(defun vector->matrix (nv nrow ncol &optional (kind :dense))
  "Copy numeric vector to a matrix of matching size."
  (let ((elements (elements nv)))
    (assert (= (length elements) (* nrow ncol)))
    (make-matrix* (lla-type nv) nrow ncol elements
                  :kind kind :shared-p (shared-p nv))))

(defun vector->column (nv)
  "Return vector as a nx1 dense matrix of the same type."
  (let ((elements (elements nv)))
    (make-matrix* (lla-type nv) (length elements) 1 elements
                  :shared-p (shared-p nv))))

(defun vector->row (nv)
  "Return vector as a 1xn dense matrix of the same type."
  (let ((elements (elements nv)))
    (make-matrix* (lla-type nv) 1 (length elements) elements
                  :shared-p (shared-p nv))))


(defun xsimilar% (rank object)
  ;; not exported
  (let ((lla-type (lla-type object))
        (rank (if (eq t rank) (xrank object) rank)))
    (ecase rank
      (1 (list (nv-class lla-type)))
      (2 `(dense-matrix :lla-type ,lla-type)))))

(defmethod xsimilar (rank (nv numeric-vector))
  (xsimilar% rank nv))

(defmethod xsimilar (rank (m dense-matrix-like))
  (xsimilar% rank m))

;; (defmethod take ((class (eql 'dense-matrix)) object &key force-copy-p options)
;;   (declare (ignore force-copy-p))
;;   (bind (((&key (lla-type :double)) options)
;;          (dims (xdims object)))
;;     (unless (= (length dims) 2)
;;       (error "OBJECT needs to have 2 dimensions for conversion into matrix."))
;;     (bind (((nrow ncol) dims)
;;            (contents (make-array (* nrow ncol)
;;                                  :element-type (lla-type->lisp-type lla-type)))
;;            (index 0))
;;       (dotimes (col ncol)
;;         (dotimes (row nrow)
;;           (setf (aref contents index) (xref object row col))
;;           (incf index)))
;;       (make-matrix nrow ncol lla-type :initial-contents contents
;;                    :type 'dense-matrix))))

(defmethod xcreate ((class (eql 'dense-matrix)) dimensions &optional options)
  (unless (= (length dimensions) 2)
    (error "Exactly 2 dimensions are needed for a matrix."))
  (bind (((nrow ncol) dimensions)
         ((&key (lla-type :double)) options))
    (make-matrix  nrow ncol lla-type)))



;; (defun take-dense-matrix-like (class matrix force-copy-p lla-type)
;;   "INTERNAL function.  Use data in a dense-matrix matrix in another
;;   dense-matrix-like object (with given class), copying and converting
;;   if necessary."
;;   (with-slots (nrow ncol data) matrix
;;     (make-instance class :nrow nrow :ncol ncol
;;                    :data (if (or force-copy-p (not (eq lla-type (lla-type matrix))))
;;                              (nv-copy-convert data lla-type)
;;                              (nv-copy data)))))

;; (defmacro define-dense-matrix-like-take (class)
;;   "Define a take method for class from dense-matrices."
;;   `(defmethod take ((class (eql ',class)) (matrix dense-matrix-like)
;;                     &key force-copy-p options)
;;      ;; set zeros
;;      (set-restricted matrix)
;;      ;; create new structure
;;      (bind (((&key (lla-type (lla-type matrix))) options))
;;        (take-dense-matrix-like ',class matrix force-copy-p lla-type))))

;; (define-dense-matrix-like-take dense-matrix)

;;;; transpose
;;;
;;; Methods should take care of returning the correct result type (eg
;;; lower-triangular into upper-triangular, etc).  Helper function
;;; transpose% can be used for implementation.

(defgeneric transpose% (matrix transposed-matrix-kind)
  (:documentation "Return the transpose of MATRIX, which will be of
kind TRANSPOSED-MATRIX-KIND.  RESTRICTED-SET-P is propagated if
meaningful, but SET-RESTRICTED is *NOT* enforced (so the caller has to
decide whether to enforce it).  Meant to be used as a helper function,
*NOT EXPORTED*."))
(expand-for-lla-types (lla-type)
  `(defmethod transpose% ((matrix ,(nv-class lla-type)) transposed-matrix-kind)
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
                 (aref elements elements-i))))
       (make-matrix* ,lla-type ncol nrow transposed :kind transposed-matrix-kind
                     :restricted-set-p (restricted-set-p matrix)))))

(defgeneric transpose (matrix)
  (:documentation "Return the transpose of a matrix.")
  (:method ((matrix dense-matrix))
    (transpose% matrix :dense))
  (:method ((matrix upper-triangular-matrix))
    (transpose% matrix :lower-triangular))
  (:method ((matrix lower-triangular-matrix))
    (transpose% matrix :upper-triangular))
  (:method ((matrix hermitian-matrix))
    (set-restricted matrix)
    (case (lla-type matrix)
      ((:complex-single :complex-double)
         (transpose% matrix 'hermitian-matrix))
      (otherwise
         (copy-matrix matrix)))))


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

(defmethod factorization-component ((mf qr-factorization) (component (eql :R)))
  (bind (((:slots-read-only nrow ncol elements) mf))
    (matrix-from-first-rows elements (lla-type mf) nrow ncol ncol :upper-triangular)))

(defmethod factorization-component ((mf cholesky-factorization) (component (eql :R)))
  (copy-matrix mf :kind :lower-triangular))


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

     
