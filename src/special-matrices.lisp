;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package :lla)

;;;; base structure and representation-independent functions

(defstruct (wrapped-matrix (:constructor make-wrapped-matrix (elements)))
  "A matrix that has some special structure (eg triangular,
symmetric/hermitian).  ELEMENTS is always a matrix."
  (elements nil :type matrix))

(defmethod elements ((matrix wrapped-matrix))
  (wrapped-matrix-elements matrix))

(defmethod diagonal ((matrix wrapped-matrix) &key copy?)
  (diagonal (elements matrix) :copy? copy?))

(defmethod emap-dimensions ((wrapped-matrix wrapped-matrix))
  (array-dimensions (elements wrapped-matrix)))

(defmethod emap-next ((wrapped-matrix wrapped-matrix) dimensions)
  (emap-next (as-array wrapped-matrix) dimensions))

(defmethod stack-into ((wrapped-matrix wrapped-matrix)
                       h? result cumulative-index)
  (stack-into (as-array wrapped-matrix) h? result cumulative-index))

;;;; mean accumulator

(defstruct (wrapped-matrix-mean-accumulator
             (:constructor wrapped-matrix-mean-accumulator% (mean kind))
             (:include array-mean-accumulator))
    "Accumulator for wrapped matrices.  Save type as extra information."
  (kind nil))

(defun wrapped-matrix-mean-accumulator (wrapped-matrix)
  "Create a conforming wrapped-matrix-mean-accumulator."
  (wrapped-matrix-mean-accumulator%
   (make-array (array-dimensions (elements wrapped-matrix))
               :initial-element 0d0)
   (matrix-kind wrapped-matrix)))

(define-conforming-accumulator (mean (matrix wrapped-matrix))
  (wrapped-matrix-mean-accumulator matrix))

(defmethod add :after ((accumulator wrapped-matrix-mean-accumulator) object)
  (let+ (((&structure wrapped-matrix-mean-accumulator- kind) accumulator))
    (when (and kind (not (eq kind (matrix-kind object))))
      (setf kind nil))))

(defmethod mean ((accumulator wrapped-matrix-mean-accumulator))
  (let+ (((&structure wrapped-matrix-mean-accumulator- mean kind) accumulator))
    (if kind
        (convert-matrix kind mean)
        mean)))

;;;; matrix creation and element access

(defgeneric matrix-kind (matrix)
  (:method ((matrix array))
    (check-type matrix matrix)
    'dense))

(declaim (inline make-matrix%))
(defun make-matrix% (nrow ncol element-type initial-element)
  "Internal function used by make-matrix methods."
  (make-array (list nrow ncol)
              :element-type element-type
              :initial-element initial-element))

(defgeneric make-matrix (kind nrow ncol
                         &key element-type initial-element)
  (:documentation "Create matrix of ")
  (:method ((kind (eql 'dense)) nrow ncol
            &key (element-type t) (initial-element (coerce 0 element-type)))
    (make-matrix% nrow ncol element-type initial-element)))

(defgeneric convert-matrix (kind object &key copy?)
  (:documentation "Convert object to a matrix of the given KIND.  May share
  structure unless COPY?.")
  (:method ((kind (eql 'dense)) object &key copy?)
    (as-array object :copy? copy?)))

(defgeneric mref (matrix row col)
  (:documentation "Element accessor for matrices.  When second value is true,
  element is constant or calculated from other elements.")
  (:method ((array array) row col)
    (aref array row col)))

(defgeneric (setf mref) (value matrix row col)
  (:documentation "Element accessor for matrices.")
  (:method (value (array array) row col)
    (setf (aref array row col) value)))

(define-condition mref-setting-readonly (error)
  ()
  (:documentation "Trying to set a constant element."))

(defgeneric represented-element? (kind row col)
  (:documentation "Return a boolean, indicating whether the element is
  represented directly.  Only defined for subclasses of WRAPPED-MATRIX."))

;;;; special matrix types

(defmacro define-special-matrix (type documentation designator
                                 represented-element? non-represented-element
                                 masked-element-string
                                 &key elements-check)
  "Define a special matrix class and methods for represents elements in
wrapped-matrix object, using a 2d array.  ROW, COL and ELEMENTS are bound to
row/column indexes and the 2d array for REPRESENTED-ELEMENT? (should evaluate
to boolean) and NON-REPRESENTED-ELEMENT (should evaluate to the value).
MASKED-ELEMENT-STRING should contain the string used for printing
non-represented elements."
  (let* ((make (symbolicate '#:make- type))
         (make-internal (if elements-check
                            (symbolicate make #\%)
                            make)))
    `(progn
       (defstruct (,type (:include wrapped-matrix)
                         (:constructor ,make-internal (elements)))
         ,documentation)
       ,@(when elements-check
           `((defun ,make (elements)
               ,elements-check
               (,make-internal elements))))
       (defmethod make-matrix ((kind (eql ',designator)) nrow ncol
                               &key (element-type t)
                                    (initial-element (coerce 0 element-type)))
         (,make (make-matrix% nrow ncol element-type initial-element)))
       (defmethod matrix-kind ((matrix ,type))
         ',designator)
       (defmethod represented-element? ((kind (eql ',designator)) row col)
         ,represented-element?)
       (defmethod convert-matrix ((kind (eql ',designator)) object &key copy?)
         (,make (as-array object :copy? copy?)))
       (defmethod nrow ((matrix ,type))
         (array-dimension (elements matrix) 0))
       (defmethod ncol ((matrix ,type))
         (array-dimension (elements matrix) 1))
       (defmethod mref ((matrix ,type) row col)
         (let+ (((&slots-r/o elements) matrix))
           (if ,represented-element?
               (aref (elements matrix) row col)
               (values ,non-represented-element t))))
       ;; !! setf mref
       (defmethod print-object ((matrix ,type) stream)
         (print-unreadable-object (matrix stream :type t)
           (terpri stream)
           (print-matrix matrix stream ,masked-element-string)))
       (defmethod as-array ((matrix ,type) &key copy?)
         (let+ (((&slots-r/o elements) matrix))
           (if copy?
               (aprog1 (make-similar-array elements)
                 (row-major-loop ((array-dimensions elements) index row col)
                   (setf (row-major-aref it index)
                         (if ,represented-element?
                             (row-major-aref elements index)
                             ,non-represented-element))))
               (prog1 elements
                 (row-major-loop ((array-dimensions elements) index row col)
                   (unless ,represented-element?
                     (setf (row-major-aref elements index)
                           ,non-represented-element))))))))))

;;;; Triangular matrices

(defun zero-like (array)
  "Return 0 coerced to the element type of array."
  (coerce 0 (array-element-type array)))

(define-special-matrix lower-triangular-matrix
    "Lower triangular matrix.  Elements in the upper triangle are zero."
  lower
  (>= row col)
  (zero-like elements)
  ".")

(define-special-matrix upper-triangular-matrix
    "Upper triangular matrix.  Elements in the upper triangle are zero."
  upper
  (<= row col)
  (zero-like elements)
  ".")

(deftype triangular-matrix ()
  '(or lower-triangular-matrix upper-triangular-matrix))

;;;; Hermitian matrices
;;;
;;; LLA uses the class HERMITIAN-MATRIX to implement both real symmetric and
;;; complex Hermitian matrices --- as technically, real symmetric matrices are
;;; also Hermitian.  Complex symmetric matrices are NOT implemented as a
;;; special matrix type, as they don't have any special properties (eg real
;;; eigenvalues, etc).

(define-special-matrix hermitian-matrix
    "Hermitian/symmetric matrix, with elements stored in the LOWER triangle."
  hermitian
  (>= row col)
  (conjugate (aref elements col row))
  "*"
  :elements-check (assert (square? elements) () "Hermitian matrices have to be
  square."))

;;;; Convenience functions for creating matrices

(defmacro dense (element-type &body rows)
  "Macro for creating a (dense) matrix (ie a rank 2 array).  ROWS should be a
list of lists (or atoms, which are treated as lists), elements are evaluated."
  (let+ ((rows (map 'vector #'ensure-list rows))
         (nrow (length rows))
         (ncol (common rows :key #'length :error
                       "Rows don't have the same number of elements."))
         ((&once-only element-type)))
    `(make-array (list ,nrow ,ncol)
                 :element-type ,element-type
                 :initial-contents
                 (list 
                  ,@(loop for row across rows collect
                          `(list
                            ,@(loop for element in row collect
                                    `(coerce ,element ,element-type))))))))

(defun pad-left-expansion% (rows ncol)
  "Pad ragged-right rows.  Used internally to implement ragged right matrix
specifications."
  (loop for row in rows
        for length = (length row)
        for row-index from 0
        collect (aprog1 (make-sequence 'list ncol :initial-element 0)
                  (replace it row :start1 0 :end1 (1+ row-index)))))

(defmacro lower (element-type &body rows)
  "Macro for creating a lower triangular matrix.  ROWS should be a list of
lists, elements are evaluated.  Masked elements (above the diagonal) are
ignored at the expansion, rows which don't have enough elements are padded
with zeros."
  `(make-lower-triangular-matrix
    (dense ,element-type
      ,@(pad-left-expansion% (mapcar #'ensure-list rows)
                             (reduce #'max rows :key #'length)))))

(defmacro hermitian (element-type &body rows)
  "Macro for creating a lower triangular matrix.  ROWS should be a list of
lists, elements are evaluated.  Masked elements (above the diagonal) are
ignored at the expansion, rows which don't have enough elements are padded
with zeros."
  `(make-hermitian-matrix
    (dense ,element-type
      ,@(pad-left-expansion% (mapcar #'ensure-list rows)
                             (max (length rows)
                                  (reduce #'max rows :key #'length))))))

(defmacro upper (element-type &body rows)
  "Macro for creating an upper triangular matrix.  ROWS should be a list of
lists, elements are evaluated.  Masked elements (below the diagonal) are
ignored at the expansion."
  `(make-upper-triangular-matrix
    (dense ,element-type
      ,@(loop for row-index from 0
              for row in rows
              collect (loop for column-index from 0
                            for element in (ensure-list row)
                            collect (if (< column-index row-index)
                                        0
                                        element))))))

(defun vec (element-type &rest elements)
  "Return a vector with elements coerced to ELEMENT-TYPE."
  (map `(simple-array ,element-type (*))
       (lambda (element) (coerce element element-type))
       elements))

;;;; Diagonal matrices

(defstruct (diagonal (:constructor make-diagonal (elements)))
  "Diagonal matrix.  The elements in the diagonal are stored in a vector."
  (elements nil :type vector))

(defmethod elements ((diagonal diagonal))
  (diagonal-elements diagonal))

(defmethod emap-dimensions ((diagonal diagonal))
  (let ((length (length (elements diagonal))))
    (list length length)))

(defmethod stack-into ((diagonal diagonal)
                       h? result cumulative-index)
  (stack-into (as-array diagonal) h? result cumulative-index))

(defmethod nrow ((diagonal diagonal))
  (length (elements diagonal)))

(defmethod ncol ((diagonal diagonal))
  (length (elements diagonal)))

(defmethod mref ((diagonal diagonal) row col)
  (if (= row col)
      (aref (elements diagonal) row)
      (values 0 t)))

(defmethod print-object ((diagonal diagonal) stream)
  (print-unreadable-object (diagonal stream :type t)
    (terpri stream)
    (print-matrix diagonal stream ".")))

(defmethod as-array ((diagonal diagonal) &key copy?)
  (declare (ignore copy?))
  (let+ (((&slots-r/o elements) diagonal)
         (n (length elements)))
    (aprog1 (make-similar-array elements :dimensions (list n n)
                                :initial-element (zero-like elements))
      (dotimes (i n) (setf (aref it i i) (aref elements i))))))

(defun diag (element-type &rest elements)
  "Return a DIAGONAL with elements coerced to ELEMENT-TYPE."
  (make-diagonal (apply #'vec element-type elements)))

;;;; elementwise operations

(defmacro define-elementwise-with-constant (type
                                            &optional (functions '(e2* e2/)))
  "Define binary elementwise operations for FUNCTION for all subclasses of
wrapped-elements."
  (let ((make (symbolicate '#:make- type)))
    `(progn
       ,@(loop :for function :in functions
               :collect
               `(defmethod ,function ((a ,type) (b number))
                  (,make (,function (elements a) b)))
               :collect
               `(defmethod ,function ((a number) (b ,type))
                  (,make (,function a (elements b))))))))
  
(defmacro define-elementwise-same-class (type
                                         &optional (functions '(e2+ e2- e2*)))
  "Define binary elementwise operations for FUNCTION for two arguments of the
same class."
  `(progn
     ,@(loop for function in functions collect
             `(defmethod ,function ((a ,type) (b ,type))
                (,(symbolicate '#:make- type)
                 (,function (elements a) (elements b)))))))

(defmacro define-elementwise-univariate 
    (type &optional (functions '(e1- e1/ eexp elog esqrt)))
  "Define unary elementwise operations for FUNCTION for all subclasses of
wrapped-elements."
  `(progn
     ,@(loop :for function :in functions
             :collect
             `(defmethod ,function ((a ,type))
                (,(symbolicate '#:make- type) (,function (elements a)))))))

(define-elementwise-with-constant lower-triangular-matrix)
(define-elementwise-with-constant upper-triangular-matrix)
(define-elementwise-with-constant hermitian-matrix)
(define-elementwise-with-constant diagonal)

(define-elementwise-same-class lower-triangular-matrix)
(define-elementwise-same-class upper-triangular-matrix)
(define-elementwise-same-class hermitian-matrix)
(define-elementwise-same-class diagonal)

(define-elementwise-univariate lower-triangular-matrix)
(define-elementwise-univariate upper-triangular-matrix)
(define-elementwise-univariate hermitian-matrix)
(define-elementwise-univariate diagonal)

;;;; transpose

(defmethod transpose ((matrix lower-triangular-matrix) &key copy?)
  (make-upper-triangular-matrix (transpose (elements matrix) :copy? copy?)))
(defmethod transpose* ((matrix lower-triangular-matrix) &key copy?)
  (make-upper-triangular-matrix (transpose* (elements matrix) :copy? copy?)))

(defmethod transpose ((matrix upper-triangular-matrix) &key copy?)
  (make-lower-triangular-matrix (transpose (elements matrix) :copy? copy?)))
(defmethod transpose* ((matrix upper-triangular-matrix) &key copy?)
  (make-lower-triangular-matrix (transpose* (elements matrix) :copy? copy?)))

(defmethod transpose ((matrix hermitian-matrix) &key copy?)
  (make-hermitian-matrix (transpose (as-array matrix) :copy? copy?)))
(defmethod transpose* ((matrix hermitian-matrix) &key copy?)
  (if copy?
      (make-hermitian-matrix (copy-array (elements matrix)))
      matrix))

(defmethod transpose ((diagonal diagonal) &key copy?)
  (if copy?
      (make-diagonal (copy-array (elements diagonal)))
      diagonal))
(defmethod transpose* ((diagonal diagonal) &key copy?)
  (declare (ignore copy?))
  (make-diagonal (econjugate (elements diagonal))))

;;;; sub
;;; 
;;; (sub wrapped-matrix index-specification nil) returns a matrix of the same
;;; class (or a scalar), using index-specification along both dimensions

(defmethod sub ((matrix wrapped-matrix) &rest index-specifications)
  (destructuring-bind (is0 &optional is1) index-specifications
    (if is1
        (sub (as-array matrix) is0 is1)
        (let ((submatrix (sub (elements matrix) is0 is0)))
          (if (arrayp submatrix)
              (convert-matrix (matrix-kind matrix) submatrix)
              submatrix)))))

;;; ==

(defmethod == ((a wrapped-matrix) (b wrapped-matrix)
               &optional (tolerance *==-tolerance*))
    (and (equal (matrix-kind a) (matrix-kind b))
         (== (as-array a) (as-array b) tolerance)))

(defmethod == ((a diagonal) (b diagonal)
               &optional (tolerance *==-tolerance*))
    (== (elements a) (elements b) tolerance))

(defgeneric as-matrix (object)
  (:documentation "Return OBJECT as an ARRAY if it has rank 2, a
  wrapped-matrix or a diagonal.")
  (:method ((array array))
    (assert (= 2 (array-rank array)))
    array)
  (:method ((matrix wrapped-matrix))
    matrix)
  (:method ((diagonal diagonal))
    diagonal))
