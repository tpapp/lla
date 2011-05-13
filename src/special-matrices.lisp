(in-package :lla)

(defclass wrapped-matrix ()
  ((elements :accessor elements :initarg :elements :type matrix
             :documentation "The representation depends on the class, use
           as-array to convert to a regular matrix."))
  (:documentation "A matrix that has some special structure (eg triangular,
  symmetric/hermitian).  Elements are always a matrix."))

(defmethod initialize-instance :after ((wrapped-matrix wrapped-matrix) 
                                       &key &allow-other-keys)
  (check-type (elements wrapped-matrix) matrix))

(defun make-array-using-contents% (dimensions initial-contents element-type
                                   copy?)
  (typecase initial-contents
    (null
       (make-array* dimensions (aif element-type it t) initial-contents))
    (number
       (make-array* dimensions (aif element-type it t) initial-contents))
    (array
       (when dimensions
         (assert (equal (ensure-list dimensions)
                        (array-dimensions initial-contents))
                 () "Dimension mismatch."))
       (if element-type
           (maybe-convert-array* initial-contents element-type copy?)
           (clnu::maybe-copy-array initial-contents copy?)))
    (otherwise (make-array-using-contents% dimensions 
                                           (as-array initial-contents)
                                           element-type copy?))))

(defmethod emap-dimensions ((wrapped-matrix wrapped-matrix))
  (array-dimensions (elements wrapped-matrix)))

(defmethod stack-into ((wrapped-matrix wrapped-matrix)
                       h? result cumulative-index)
  (stack-into (as-array wrapped-matrix) h? result cumulative-index))

(defmethod mean-accumulator ((first-element wrapped-matrix) sequence)
  (let ((type (type-of first-element))
        (array-accumulator (mean-accumulator (elements first-element)
                                             sequence)))
    (lambda (&optional x)
      (if x
          (progn
            (when (and type (not (equal type (type-of x))))
              (setf type nil))
            (funcall array-accumulator x))
          (let ((mean (funcall array-accumulator)))
            (if type
                (make-instance type :elements mean)
                mean))))))

(defgeneric make-matrix (kind dimensions 
                              &key initial-contents element-type copy?)
  (:documentation "Create a matrix of given KIND, with DIMENSIONS (can be NIL
  if inferred from initial contents, usually another matrix-like object).")
  (:method ((kind (eql :dense)) dimensions
            &key initial-contents element-type copy?)
    (make-array-using-contents% dimensions initial-contents element-type
                                copy?)))

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

(defmacro define-special-matrix (class documentation synonyms
                                 represented-element? non-represented-element
                                 masked-element-string &key bindings)
  "Define a special matrix class and methods for represents elements in
 wrapped-matrix object, using a 2d array.  ROW, COL and ELEMENTS are bound to
 row/column indexes and the 2d array for REPRESENTED-ELEMENT? (should evaluate
 to boolean) and NON-REPRESENTED-ELEMENT (should evaluate to the value).
 MASKED-ELEMENT-STRING should contain the string used for printing
 non-represented elements."
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
                (make-matrix ',class dimensions :copy? copy?
                             :element-type element-type
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
                         ,non-represented-element)))))))))

;;; Triangular matrices

(define-special-matrix lower-triangular-matrix
    "Lower triangular matrix.  Elements in the upper triangle are zero."
  (:lower :lower-triangular)
  (>= row col)
  zero
  "."
  :bindings ((zero (zero* elements))))

(define-special-matrix upper-triangular-matrix
    "Upper triangular matrix.  Elements in the upper triangle are zero."
  (:upper :upper-triangular)
  (<= row col)
  zero
  "."
  :bindings ((zero (zero* elements))))

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

(declaim (inline hermitian-orientation))
(defun hermitian-orientation (library)
  "Return a constant that denotes the orientation of hermitian matrices."
  (ecase library
    (:blas :CBLASLOWER)
    (:lapack +l+)))

(defmethod initialize-instance :after ((hermitian-matrix hermitian-matrix)
                                       &key &allow-other-keys)
  (check-type (elements hermitian-matrix) (and matrix (satisfies square?))))

;;; Diagonal matrices

(defclass diagonal-matrix ()
  ((elements :accessor elements :initarg :elements :type vector
             :documentation "The representation depends on the class, use
           as-array to convert to a regular matrix."))
  (:documentation
   "Diagonal matrix.  The elements in the diagonal are stored in a vector."))

(defmethod initialize-instance :after ((diagonal-matrix diagonal-matrix)
                                       &key &allow-other-keys)
  (check-type (elements diagonal-matrix) vector))

(defmethod emap-dimensions ((diagonal-matrix diagonal-matrix))
  (let ((length (length (elements diagonal-matrix))))
    (list length length)))

(defmethod stack-into ((diagonal-matrix diagonal-matrix)
                       h? result cumulative-index)
  (stack-into (as-array diagonal-matrix) h? result cumulative-index))

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
                                :initial-element (zero* elements))
      (dotimes (i n) (setf (aref it i i) (aref elements i))))))

;;; elementwise operations

(defmacro define-elementwise-with-constant (class &optional (functions '(e2* e2/)))
  "Define elementwise operations for FUNCTION for all subclasses of
wrapped-elements.  "
  `(progn
     ,@(loop :for function :in functions
             :collect
             `(defmethod ,function ((a ,class) (b number))
                (make-instance ',class :elements (,function (elements a) b)))
             :collect
             `(defmethod ,function ((a number) (b ,class))
                (make-instance ',class :elements (,function a (elements b)))))))
  
(defmacro define-elementwise-same-class (class &optional (functions '(e2+ e2- e2*)))
  `(progn
     ,@(loop for function in functions collect
             `(defmethod ,function ((a ,class) (b ,class))
                (make-instance ',class
                               :elements (,function (elements a) (elements b)))))))

(define-elementwise-with-constant lower-triangular-matrix)
(define-elementwise-with-constant upper-triangular-matrix)
(define-elementwise-with-constant hermitian-matrix)
(define-elementwise-with-constant diagonal-matrix)

(define-elementwise-same-class lower-triangular-matrix)
(define-elementwise-same-class upper-triangular-matrix)
(define-elementwise-same-class hermitian-matrix)
(define-elementwise-same-class diagonal-matrix)

;;; transpose

(defmethod transpose ((matrix lower-triangular-matrix) &key copy?)
  (make-instance 'upper-triangular-matrix
                 :elements (transpose (elements matrix) :copy? copy?)))
(defmethod transpose* ((matrix lower-triangular-matrix) &key copy?)
  (make-instance 'upper-triangular-matrix
                 :elements (transpose* (elements matrix) :copy? copy?)))

(defmethod transpose ((matrix upper-triangular-matrix) &key copy?)
  (make-instance 'lower-triangular-matrix
                 :elements (transpose (elements matrix) :copy? copy?)))
(defmethod transpose* ((matrix upper-triangular-matrix) &key copy?)
  (make-instance 'lower-triangular-matrix
                 :elements (transpose* (elements matrix) :copy? copy?)))

(defmethod transpose ((matrix hermitian-matrix) &key copy?)
  (make-instance 'hermitian-matrix
                 :elements (transpose (as-array matrix) :copy? copy?)))
(defmethod transpose* ((matrix hermitian-matrix) &key copy?)
  (if copy?
      (make-instance 'hermitian-matrix :elements (copy-array (elements matrix)))
      matrix))

(defmethod transpose ((diagonal diagonal-matrix) &key copy?)
  (if copy?
      (make-instance 'diagonal-matrix :elements (copy-array (elements diagonal)))
      diagonal))
(defmethod transpose* ((diagonal diagonal-matrix) &key copy?)
  (declare (ignore copy?))
  (make-instance 'diagonal-matrix :elements (econjugate (elements diagonal))))

;;; sub
;;; 
;;; (sub wrapped-matrix index-specification nil) returns a matrix of the same class
;;; (or a scalar), using index-specification along both dimensions

(defmethod sub ((matrix wrapped-matrix) &rest index-specifications)
  (destructuring-bind (is0 &optional is1) index-specifications
    (if is1
        (sub (as-array matrix) is0 is1)
        (let ((submatrix (sub (elements matrix) is0 is0)))
          (if (arrayp submatrix)
              (make-instance (class-of matrix) :elements submatrix)
              submatrix)))))

;;; ==

(defmethod == ((a wrapped-matrix) (b wrapped-matrix)
               &optional (tolerance *==-tolerance*))
    (and (equal (type-of a) (type-of b))
         (== (as-array a) (as-array b) tolerance)))

(defmethod == ((a diagonal-matrix) (b diagonal-matrix)
               &optional (tolerance *==-tolerance*))
    (== (elements a) (elements b) tolerance))
