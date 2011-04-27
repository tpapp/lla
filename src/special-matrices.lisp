(in-package :lla)

(defclass wrapped-matrix ()
  ((elements :accessor elements :initarg :elements :type elements
           :documentation "The representation depends on the class, use as-array to
           convert to a regular matrix."))
  (:documentation "A matrix that has some special structure (eg triangular,
  symmetric/hermitian)."))

(defmethod initialize-instance :after ((wrapped-matrix wrapped-matrix) 
                                       &key &allow-other-keys)
  (check-type (elements wrapped-matrix) matrix))

(defun make-array-using-contents% (dimensions initial-contents element-type copy?)
  (typecase initial-contents
    (null
     (make-array* dimensions (aif element-type it t) initial-contents))
    (number
     (make-array* dimensions (aif element-type it t) initial-contents))
    (array
     (when dimensions
       (assert (equal (ensure-list dimensions) (array-dimensions initial-contents))
               () "Dimension mismatch."))
     (if element-type
         (maybe-convert-array* initial-contents element-type copy?)
         (clnu::maybe-copy-array initial-contents copy?)))
    (otherwise (make-array-using-contents% dimensions (as-array initial-contents)
                                   element-type copy?))))

(defgeneric make-matrix (kind dimensions &key initial-contents element-type copy?)
  (:documentation "Create a matrix of given KIND, with DIMENSIONS (can be NIL if
  inferred from initial contents, usually another matrix-like object).")
  (:method ((kind (eql :dense)) dimensions &key initial-contents element-type copy?)
    (make-array-using-contents% dimensions initial-contents element-type copy?)))

(defgeneric mref (matrix row col)
  (:documentation "Element accessor for matrices.  When second value is true, element
  is constant or calculated from other elements.")
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
  (:documentation "Return a boolean, indicating whether the element is represented
  directly.  Only defined for subclasses of WRAPPED-MATRIX."))

(defmacro define-special-matrix (class documentation synonyms represented-element?
                                 non-represented-element masked-element-string &key bindings)
  "Define a special matrix class and methods for represents elements in wrapped-matrix object,
using a 2d array.  ROW, COL and ELEMENTS are bound to row/column indexes and the 2d
array for REPRESENTED-ELEMENT? (should evaluate to boolean) and
NON-REPRESENTED-ELEMENT (should evaluate to the value).  MASKED-ELEMENT-STRING should
contain the string used for printing non-represented elements."
  `(progn
     (defclass ,class (wrapped-matrix) () (:documentation ,documentation))
     (defmethod make-matrix ((kind (eql ',class)) dimensions &key
                             initial-contents element-type copy?)
       (make-instance ',class :elements
                      (make-array-using-contents% dimensions initial-contents
                                                  element-type copy?)))
     ,@(loop for synonym :in synonyms collect
             `(defmethod make-matrix ((kind (eql ',synonym)) dimensions &key
                                      initial-contents element-type copy?)
                (make-matrix ',class dimensions :copy? copy? :element-type element-type
                             :initial-contents initial-contents)))
     (defmethod represented-element? ((kind (eql ',class)) row col)
       ,represented-element?)
     ,@(loop for synonym :in synonyms :collect
             `(defmethod represented-element? ((kind (eql ',synonym)) row col)
                (represented-element? ',class row col)))
     (defmethod nrow ((matrix ,class)) (array-dimension (elements matrix) 0))
     (defmethod ncol ((matrix ,class)) (array-dimension (elements matrix) 1))
     (defmethod mref ((matrix ,class) row col)
       (bind (((:slots-r/o elements) matrix)
              ,@bindings)
         (if ,represented-element?
             (aref (elements matrix) row col)
             (values ,non-represented-element t))))
     ;; !! setf mref
     (defmethod print-object ((matrix ,class) stream)
       (print-unreadable-object (matrix stream :type t)
         (terpri stream)
         (print-matrix matrix stream ,masked-element-string)))
     (defmethod as-array ((matrix ,class) &key copy?)
       (bind (((:slots-r/o elements) matrix)
              ,@bindings)
         (if copy?
             (aprog1 (make-similar-array matrix)
               (row-major-loop (elements index row col)
                 (setf (row-major-aref it index)
                       (if ,represented-element?
                           (row-major-aref elements index)
                           ,non-represented-element))))
             (prog1 elements
               (row-major-loop (elements index row col)
                 (unless ,represented-element?
                   (setf (row-major-aref elements index)
                         ,non-represented-element)))))))))

;;; Triangular matrices

(define-special-matrix lower-triangular-matrix
    "Lower triangular matrix.  Elements in the upper triangle are zero."
  (:lower :lower-triangular)
  (>= row col)
  zero
  "."
  :bindings ((zero (zero* (array-element-type elements)))))

(define-special-matrix upper-triangular-matrix
    "Upper triangular matrix.  Elements in the upper triangle are zero."
  (:upper :upper-triangular)
  (<= row col)
  zero
  "."
  :bindings ((zero (zero* (array-element-type elements)))))

(deftype triangular-matrix ()
  '(or lower-triangular-matrix upper-triangular-matrix))

;;; Hermitian matrices
;;;
;;; LLA uses the class HERMITIAN-MATRIX to implement both real symmetric and complex
;;; Hermitian matrices --- as technically, real symmetric matrices are also
;;; Hermitian.  Complex symmetric matrices are NOT implemented as a special matrix
;;; type, as they don't have any special properties (eg real eigenvalues, etc).

(define-special-matrix hermitian-matrix
    "Hermitian/symmetric matrix, with elements stored in the lower triangle."
  (:hermitian :hermitian-matrix)
  (>= row col)
  (conjugate (aref elements col row))
  "*")

(defmethod initialize-instance :after ((hermitian-matrix hermitian-matrix)
                                       &key &allow-other-keys)
  (check-type (elements hermitian-matrix) (and matrix (satisfies square?))))

;;; Diagonal matrices

(defclass diagonal-matrix ()
  ((elements :accessor elements :initarg :elements))
  (:documentation "Diagonal matrix.  The elements in the diagonal are stored in a
  vector."))

(defmethod initialize-instance :after ((diagonal-matrix diagonal-matrix)
                                       &key &allow-other-keys)
  (check-type (elements diagonal-matrix) vector))

(defmethod make-matrix ((kind (eql 'diagonal-matrix)) dimensions 
                        &key initial-contents element-type copy?)
  (make-instance 'diagonal-matrix :elements 
                 (make-array-using-contents% dimensions initial-contents
                                             element-type copy?)))

(defmethod make-matrix ((kind (eql :diagonal)) dimensions 
                        &key initial-contents element-type copy?)
  (make-matrix 'diagonal-matrix dimensions :copy? copy? :element-type element-type
               :initial-contents initial-contents))

(defmethod nrow ((diagonal-matrix diagonal-matrix))
  (length (elements diagonal-matrix)))

(defmethod ncol ((diagonal-matrix diagonal-matrix))
  (length (elements diagonal-matrix)))

(defmethod mref ((diagonal-matrix diagonal-matrix) row col)
  (if (= row col)
      (aref (elements diagonal-matrix) row)
      (values 0 t)))

(defmethod print-object ((diagonal-matrix diagonal-matrix) stream)
  (print-unreadable-object (diagonal-matrix stream :type t)
    (terpri stream)
    (print-matrix diagonal-matrix stream ".")))

(defmethod as-array ((diagonal-matrix diagonal-matrix) &key copy?)
  (declare (ignore copy?))
  (bind (((:slots-r/o elements) diagonal-matrix)
         (n (length elements)))
    (aprog1 (make-similar-array elements :dimensions (list n n)
                                :initial-element (zero* (array-element-type elements)))
      (dotimes (i n) (setf (aref it i i) (aref elements i))))))
