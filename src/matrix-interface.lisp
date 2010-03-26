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

(defmethod set-restricted ((matrix upper-triangular-matrix))
  (declare (optimize speed))
  (bind (((:slots-read-only nrow ncol elements) matrix)
         (zero (coerce* 0 (lla-type matrix))))
         ;; set the lower triangle (below diagonal) to 0
    (declare (fixnum nrow ncol)
             (type simple-array elements))
    (dotimes (col ncol)
      (declare (fixnum col ncol))
      (iter
        (declare (iterate:declare-variables))
        (for (the fixnum index)
             :from (1+ (cm-index2 nrow col col))
             :below (cm-index2 nrow nrow col))
        (setf (aref elements index) zero))))
  matrix)

(defmethod set-restricted ((matrix lower-triangular-matrix))
   (bind (((:slots-read-only nrow ncol elements) matrix)
         (zero (zero* (lla-type matrix))))
    ;; set the upper triangle (above diagonal) to 0
    (dotimes (col ncol)
      (iter
        (for index
          :from (cm-index2 nrow 0 col)
          :below (cm-index2 nrow col col))
        (setf (aref elements index) zero))))
  matrix)

(defmethod set-restricted ((matrix hermitian-matrix))
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
;;; Notes: XELTTYPE is already defined for the NUMERIC-VECTOR-LIKE
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
  (bind (((row col) subscripts)
         ((:accessors-read-only leading-dimension nrow ncol elements) matrix))
    (check-index row nrow)
    (check-index col ncol)
    (aref elements (cm-index2 leading-dimension row col))))

(defun matrix-setf-xref% (value matrix subscripts)
  "Setter function for matrices.  Meant to be inlined.  Not exported."
  (bind (((row col) subscripts)
         ((:accessors-read-only leading-dimension nrow ncol elements) matrix))
    (check-index row nrow)
    (check-index col ncol)
    (setf (aref elements (cm-index2 leading-dimension row col))
          value)))

(defmethod xref ((matrix dense-matrix-like) &rest subscripts)
  (matrix-xref% matrix subscripts))

(defmethod (setf xref) (value (matrix dense-matrix-like) &rest subscripts)
  (matrix-setf-xref% value matrix subscripts))

;;; xref for upper triangular matrices

(defmethod xref ((matrix upper-triangular-matrix) &rest subscripts)
  (bind (((row col) subscripts))
    (with-slots (nrow ncol elements) matrix
      (check-index row nrow)
      (check-index col ncol)
      (if (<= row col)
          (aref elements (cm-index2 nrow row col))
          (zero* (lla-type matrix))))))

(defmethod (setf xref) (value (matrix upper-triangular-matrix) &rest subscripts)
  (bind (((row col) subscripts))
    (with-slots (nrow ncol elements) matrix
      (check-index row nrow)
      (check-index col ncol)
      (if (<= row col)
          (setf (aref elements (cm-index2 nrow row col))
                value)
          (if (zerop value)
              value
              (error 'xref-setting-readonly))))))

;;; xref for lower triangular matrices

(defmethod xref ((matrix lower-triangular-matrix) &rest subscripts)
  (bind (((row col) subscripts))
    (with-slots (nrow ncol elements) matrix
      (check-index row nrow)
      (check-index col ncol)
      (if (>= row col)
          (aref elements (cm-index2 nrow row col))
          (zero* (lla-type matrix))))))

(defmethod (setf xref) (value (matrix lower-triangular-matrix) &rest subscripts)
  (bind (((row col) subscripts))
    (with-slots (nrow ncol elements) matrix
      (check-index row nrow)
      (check-index col ncol)
      (if (>= row col)
          (setf (aref elements (cm-index2 nrow row col)) value)
          (if (zerop value)
              value
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
  (bind (((row col) subscripts))
    (with-slots (nrow ncol elements) matrix
      (check-index row nrow)
      (check-index col ncol)
      (if (<= row col)
          (setf (aref elements (cm-index2 nrow row col)) value)
          (setf (aref elements (cm-index2 nrow col row)) (conjugate value))))))

;;;; matrix creation

(declaim (inline make-matrix*))
(defun make-matrix* (lla-type nrow ncol elements &key (kind :dense))
  "Create a matrix with given ELEMENTS, TYPE, LLA-TYPE and dimensions.
Note that there is no type checking, and elements are not copied: this
is effectively shorthand for a MAKE-INSTANCE call.  For internal use,
not exported."
  (make-instance (matrix-class kind) :nrow nrow :ncol ncol
                 :lla-type lla-type :elements elements))

(defun make-matrix (lla-type nrow ncol &key (kind :dense) (initial-element 0))
  "Create a matrix with given parameters, optionally initialized with
INITIAL-ELEMENTs."
  (make-matrix* lla-type nrow ncol
                (make-nv-elements lla-type (* nrow ncol) initial-element)
                :kind kind))

(defun matrix-elements-from-sequence (ncol lla-type sequence)
  "Return a lisp vector comforming to LLA-TYPE which has initial
contents derived from 2d Lisp ARRAY.  Return (VALUES ELEMENTS NROW NCOL).
If NCOL is zero, return a row matrix."
  (bind ((length (length sequence))
         ((:values nrow ncol) (if (zerop ncol)
                                  (values 1 length)
                                  (bind (((:values nrow remainder) (floor length ncol)))
                                    (unless (zerop remainder)
                                      (error "Length of sequence (~A) is not ~
                                              a multiple of ncol (~A)." length ncol))
                                    (values nrow ncol))))
	 (lisp-type (lla-type->lisp-type lla-type))
	 (elements (make-nv-elements lla-type length))
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
    (values elements nrow ncol)))

(defun create-matrix (ncol initial-contents &key (kind :dense) lla-type)
  "Create matrix of given TYPE vector with given initial contents (a
sequence).  Unless LLA-TYPE is given, it is inferred from the
elements.  NCOL gives the number of columns, while the number of rows
is inferred from the length of the sequence.  An error is signalled if
there are remainder elements.  Elements corresponding to restricted
elements are just ignored.

Usage note: This is a convenience function for easily creation of
matrices.  Also see *force-float*."
  (bind ((lla-type (infer-lla-type lla-type initial-contents))
         ((:values elements nrow ncol)
          (matrix-elements-from-sequence ncol lla-type initial-contents)))
    (make-matrix* lla-type nrow ncol elements :kind kind)))

(defun copy-matrix% (matrix &key (kind (matrix-kind matrix))
                     (destination-type (lla-type matrix)) (copy-p nil))
  "Copy or convert matrix to the given kind and destination-type.
Copying is forced when COPY-P.  Important usage note: if you want a
matrix with restricted elements, it is advisable to set copy-p,
otherwise a call to SET-RESTRICTED hidden somewhere might change the
original matrix, without your intention.  Also, SET-RESTRICTED is not
called by this function.  That's why this function is not exported.
Use with caution."
  (make-matrix* (lla-type matrix) (nrow matrix) (ncol matrix)
                (copy-elements% matrix destination-type copy-p)
                :kind kind))

(defun reshape (nv nrow ncol &key (kind :dense) (copy-p nil))
  "Copy numeric vector (can be a matrix) to a matrix of matching
size."
  (let ((elements (elements nv)))
    (assert (= (length elements) (* nrow ncol)))
    (make-matrix* (lla-type nv) nrow ncol
                  (if copy-p (copy-seq elements) elements)
                  :kind kind)))

(defgeneric vector->column (vector &key copy-p)
  (:documentation "Return vector as a nx1 dense matrix of the same type."))

(defmethod vector->column ((nv numeric-vector) &key (copy-p nil))
  (let ((elements (elements nv)))
    (make-matrix* (lla-type nv) (length elements) 1
                  (if copy-p (copy-seq elements) elements))))

(defgeneric vector->row (vector &key copy-p)
  (:documentation "Return vector as a 1xn dense matrix of the same type."))

(defmethod vector->row ((nv numeric-vector) &key (copy-p nil))
  "Return vector as a 1xn dense matrix of the same type."
  (let ((elements (elements nv)))
    (make-matrix* (lla-type nv) 1 (length elements)
                  (if copy-p (copy-seq elements) elements))))

(defmethod make-load-form ((matrix dense-matrix-like) &optional environment)
  (declare (ignore environment))
  `(make-matrix* ,(lla-type matrix) ,(nrow matrix) ,(ncol matrix) 
                 ,(make-elements-load-form matrix) :kind ,(matrix-kind matrix)))
