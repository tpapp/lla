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

(expand-for-lla-types lla-type
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
  (bind (((:slots-read-only nrow ncol data) matrix)
         (vector (copy-data data)))
    ;; set the lower triangle (below diagonal) to conjugate of the
    ;; elements in the upper triangle
    (dotimes (col ncol)
      (iter
        (for row :from col :below nrow)
        (for index
          :from (cm-index2 nrow col col)
          :below (cm-index2 nrow nrow col))
        (setf (aref vector index)
              (conjugate (aref vector (cm-index2 nrow col row)))))))
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

(defmethod xref ((matrix dense-matrix) &rest subscripts)
  (bind (((row col) subscripts))
    (with-slots (nrow ncol elements) matrix
      (check-index row nrow)
      (check-index col ncol)
      (aref elements (cm-index2 nrow row col)))))

(defmethod (setf xref) (value (matrix dense-matrix) &rest subscripts)
  (ensure-unshared matrix)
  (bind (((row col) subscripts))
    (with-slots (nrow ncol elements) matrix
      (check-index row nrow)
      (check-index col ncol)
      (setf (aref elements (cm-index2 nrow row col))
	    value))))

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


;;;; matrix creation

(declaim (inline make-matrix*))
(defun make-matrix* (lla-type nrow ncol elements &key (kind 'dense)
                     (shared-p nil))
  "Create a matrix with given ELEMENTS, TYPE, LLA-TYPE and dimensions.
Note that there is no type checking, and elements are not copied: this
is effectively shorthand for a MAKE-INSTANCE call.  For internal use,
not exported."
  (make-instance (matrix-class kind lla-type) :nrow nrow :ncol ncol
                 :elements elements :shared-p shared-p))

(defun make-matrix (nrow ncol lla-type &key (kind :dense) initial-element)
  "Create a matrix with given parameters, optionally initialized with
INITIAL-ELEMENTs (where it makes sense, for restricted matrices)."
  (let ((lisp-type (lla-type->lisp-type lla-type))
        (length (* nrow ncol)))
    (make-matrix* lla-type nrow ncol
                  (if initial-element
                      (make-array length :element-type lisp-type
                                  :initial-element (coerce initial-element lisp-type))
                      (make-array length :element-type lisp-type))
                  :kind kind)))

(defun matrix-elements-from-array (lla-type array)
  "Return a NUMERIC-VECTOR of LLA-TYPE which has initial contents
derived from 2d Lisp ARRAY.  Return (VALUES ELEMENTS NROW NCOL).
Fails if array is not 2D."
  (bind (((nrow ncol) (array-dimensions array))
	 (lisp-type (lla-type->lisp-type lla-type))
	 (elements (make-array (* nrow ncol) :element-type lisp-type)))
    (dotimes (row nrow)
      (dotimes (col ncol)
	(setf (aref elements (cm-index2 nrow row col))   ; column major
	      (coerce (aref array row col) lisp-type)))) ; row major
    (values elements nrow ncol)))

(defun create-matrix (initial-contents &key (kind :dense) lla-type)
  "Create matrix of given TYPE vector with given initial contents (a
2D array).  Unless LLA-TYPE is given, it is inferred from the
elements.  This is a convenience function for easily creation of
matrices.  Also see *force-float*."
  (bind ((lla-type (infer-lla-type lla-type (flatten-array initial-contents) *force-float*))
         ((:values elements nrow ncol) (matrix-elements-from-array lla-type initial-contents)))
    (make-matrix* lla-type nrow ncol elements :kind kind)))

(defun copy-matrix (matrix &optional (kind (matrix-kind matrix)))
  "Copy matrix.  When applicable, elements will be restricted.  Usage
note: if you want to get a numeric vector, just use COPY-NV."
  (set-restricted matrix)
  (bind (((:slots-read-only nrow ncol elements shared-p) matrix)
         (result (make-matrix* (lla-type matrix) nrow ncol elements
                               :kind kind :shared-p shared-p)))
        (when (typep result 'restricted-elements)
          (setf (restricted-set-p result) t))
        result))
  
(defun convert-matrix (matrix lla-type &optional (kind :dense))
  "Convert matrix.  When applicable, elements will be restricted.
Usage note: if you want to get a numeric vector, just use CONVERT-NV."
  (set-restricted matrix)
  (let ((result (make-matrix* lla-type (nrow matrix) (ncol matrix)
                              (convert-elements matrix lla-type)
                              :kind kind)))
    (when (typep result 'restricted-elements)
      (setf (restricted-set-p result) t))
    result))

(defun vector->matrix (nv nrow ncol &optional (kind 'dense))
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


;; (defun xsimilar% (rank object)
;;   ;; not exported
;;   (let ((lla-type (lla-type object))
;;         (rank (if (eq t rank) (xrank object) rank)))
;;     (ecase rank
;;       (1 (list (numeric-vector-class lla-type)))
;;       (2 `(dense-matrix :lla-type ,lla-type)))))

;; (defmethod xsimilar (rank (nv numeric-vector))
;;   (xsimilar% rank nv))

;; (defmethod xsimilar (rank (m dense-matrix-like))
;;   (xsimilar% rank m))



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

;; (defmethod xcreate ((class (eql 'dense-matrix)) dimensions &optional options)
;;   (unless (= (length dimensions) 2)
;;     (error "Exactly 2 dimensions are needed for a matrix."))
;;   (bind (((nrow ncol) dimensions)
;;          ((&key (lla-type :double)) options))
;;     (make-matrix class nrow ncol :lla-type lla-type)))



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
;;;;
;;;; Methods should take care of returning the correct result type (eg
;;;; lower-triangular into upper-triangular, etc).  Helper function
;;;; transpose% can be used for implementation.

;; (defun transpose% (matrix transposed-matrix-class)
;;   "Return the transpose of MATRIX, which will be of class
;; TRANSPOSED-MATRIX-CLASS.  RESTRICTED-SET-P is propagated if
;; meaningful.  Meant to be used as a helper function."
;;   (declare (optimize (speed 3)))
;;   (bind (((:slots-read-only nrow ncol data) matrix)
;;          (original-vector (nv-data data))
;;          (transposed (make-nv (length original-vector) (lla-type data)))
;;          (transposed-vector (nv-data transposed)))
;;     (declare (fixnum nrow ncol))
;;     ;; !!! could do with a bit of optimization, eg methods specialized
;;     ;; !!! on lla-type, some fixnums are still not recognized --- Tamas
;;     (dotimes (col ncol)
;;       (declare (fixnum col))
;;       (iter
;;         (declare (iterate:declare-variables))
;;         (for (the fixnum original-i) :from (cm-index2 nrow 0 col) :below (cm-index2 nrow nrow col))
;;         (for (the fixnum transposed-i) :from (cm-index2 ncol col 0) :by ncol)
;;         (setf (aref transposed-vector transposed-i)
;;               (aref original-vector original-i))))
;;     (let ((transposed-matrix (make-matrix transposed-matrix-class ncol nrow
;;                                           :lla-type (lla-type matrix)
;;                                           :initial-contents transposed
;;                                           :use-directly-p t)))
;;       (if (subtypep transposed-matrix-class 'dense-matrix-restricted)
;;           (setf (restricted-set-p transposed-matrix) (restricted-set-p matrix)))
;;       transposed-matrix)))

;; (defmethod transpose ((matrix dense-matrix))
;;   (transpose% matrix 'dense-matrix))

;; (defmethod x* ((a dense-matrix-like) (b number) &key)
;;   ;;; !!!! this assumes that all restricted elements are zeros
;;   (if (typep a 'restricted-elements)
;;       (make-instance (class-of a) :nrow (nrow a) :ncol (ncol a)
;;                      :data (x* (data a) b) :restricted-set-p (restricted-set-p a))
;;       (make-instance (class-of a) :nrow (nrow a) :ncol (ncol a)
;;                      :data (x* (data a) b))))

;; (defmethod x* ((a number) (b dense-matrix-like) &key)
;;   (x* b a))


;; (defmethod transpose ((matrix upper-triangular-matrix))
;;   (transpose% matrix 'lower-triangular-matrix))

;; (defmethod transpose ((matrix lower-triangular-matrix))
;;   (transpose% matrix 'upper-triangular-matrix))

;; (defmethod transpose ((matrix hermitian-matrix))
;;   ;; !!! could just leave real matrices alone
;;   (transpose% matrix 'hermitian-matrix))


  



;; ;;;; ?? define xarray interface for LU, remember that it has 3
;; ;;;; dimensions, and it is pretty meaningless to access L and U
;; ;;;; separately.  currently, I don't see a need for this, it doesn't
;; ;;;; integrate nicely into xarray.  Let's keep this semi-opauque
;; ;;;; (those who really want it can access ipiv and data and interpret
;; ;;;; it according to LAPACK conventions). -- Tamas

;; ;;;; !! define an initialize-instance after method for consistency check?

;; ;;;; !! define methods for factorization components

;; ;;;;
;; ;;;;  QR factorization
;; ;;;;

;; (defmethod factorization-component ((mf QR) (component (eql :R)) 
;;                                     &optional force-copy-p)
;;   (declare (ignore force-copy-p))
;;   (bind (((:slots-read-only nrow ncol data) mf))
;;     (take 'upper-triangular-matrix
;;           (matrix-from-first-rows data nrow ncol ncol))))

;; ;;;; !! define an initialize-instance after method for consistency check?


;; ;;;;
;; ;;;;  Cholesky decomposition
;; ;;;;


;; (defmethod initialize-instance :after ((mf cholesky) &key &allow-other-keys)
;;   (with-slots (nrow ncol data) mf
;;     (assert (= nrow ncol))
;;     (assert (equalp (xdims data) (list (* nrow ncol)))))
;;   mf)

;; (defmethod factorization-component ((mf cholesky) component &optional force-copy-p)
;;   (declare (ignore force-copy-p))
;;   (take 'lower-triangular-matrix mf))

