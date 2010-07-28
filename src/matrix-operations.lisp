(in-package :lla)

;;;; extraction methods for factorization components and related
;;;; utility functions

(defun matrix-from-first-rows (vector nrow ncol leading-dimension
                               &optional (kind :dense))
  "Extract & return (as KIND, default is :DENSE) the first NROW rows
of an LEADING-DIMENSIONxNCOL matrix, given as a Lisp vector with
column-major indexing.  NOTE: needed to interface to LAPACK routines
like xGELS."
  (let ((result (make-similar-vector vector (* nrow ncol))))
    (copy-columns nrow ncol 
                  vector 0 leading-dimension
                  result 0 nrow)
    (make-matrix% nrow ncol result :kind kind)))

;;;; stacking

(defgeneric dimensions-as-matrix% (object horizontal?)
  (:documentation "Return dimensions when used as a matrix for
  stacking.  HORIZONTAL? determines the direction for stacking.
  Return (values NROW NCOL)")
  (:method ((object dense-matrix-like) horizontal?) 
    (values (nrow object) (ncol object)))
  (:method ((object diagonal) horizontal?)
    (let ((l (length (elements object))))
      (values l l)))
  (:method ((object vector) horizontal?)
    (let ((l (length object)))
      (if horizontal?
          (values l 1)
          (values 1 l))))
  (:method ((object array) horizontal?)
    (bind (((nrow ncol) (array-dimensions object)))
      (values nrow ncol))))

(defgeneric fill-matrix-elements-using% (object matrix-elements offset
                                                leading-dimension
                                                horizontal?)
  (:documentation "Fill MATRIX-ELEMENTS (starting at OFFSET, with
  LEADING-DIMENSION) using OBJECT.  HORIZONTAL? determines the
  direction for stacking.  Return no values.")
  (:method ((vector vector) matrix-elements offset leading-dimension
            horizontal?)
    (if horizontal?
        (copy-elements vector 0 
                       matrix-elements offset (length vector))
        (copy-columns 1 (length vector)
                      vector 0 1
                      matrix-elements offset leading-dimension)))
  (:method ((diagonal diagonal) matrix-elements offset
            leading-dimension horizontal?)
    (bind (((:slots-r/o elements) diagonal))
      ;; dirty hack: we use 1+leading-dimension to "slide" the diagonal
      (copy-columns 1 (length elements)
                    elements 0 1
                    matrix-elements offset (1+ leading-dimension))))
  (:method ((matrix dense-matrix-like) matrix-elements offset
            leading-dimension horizontal?)
    (set-restricted matrix)
    (bind (((:slots-r/o elements nrow ncol) matrix))
      (copy-columns nrow ncol
                    elements 0 nrow
                    matrix-elements offset leading-dimension)))
  (:method ((matrix array) matrix-elements offset leading-dimension
            horizontal?)
    (fill-matrix-elements-using% (as-matrix matrix) matrix-elements
                                 offset
                                 leading-dimension horizontal?)))

(defun stack (direction objects)
  "Stack objects (a list) in the given direction (:horizontal/:h
or :vertical/:v)."
  (bind ((horizontal? (ecase direction
                        ((:horizontal :h) t)
                        ((:vertical :v) nil)))
         ((:values nrows ncols)
          (iter
            (for object :in objects)
            (for (values nrow ncol) := 
                 (dimensions-as-matrix% object horizontal?))
            (collecting nrow :into nrows)
            (collecting ncol :into ncols)
            (finally
             (return (values nrows ncols))))))
    (assert (apply #'= (if horizontal? nrows ncols)) ()
            "Dimensions don't match.")
    (bind ((offset 0)
           (lla-type (reduce #'common-lla-type objects
                             :key (lambda (obj)
                                    (array-lla-type 
                                     (if (typep obj 'elements%)
                                         (elements obj)
                                         obj))))))
      (if horizontal?
          (bind ((nrow (first nrows))
                 (result (make-matrix nrow (reduce #'+ ncols) lla-type))
                 ((:slots-r/o elements) result))
            (iter
              (for ncol% :in ncols)
              (for object :in objects)
              (fill-matrix-elements-using% object elements
                                           offset nrow t)
              (incf offset (* ncol% nrow)))
            result)
          (bind ((nrow (reduce #'+ nrows))
                 (ncol (first ncols))
                 (result (make-matrix nrow ncol lla-type))
                 ((:slots-r/o elements) result))
            (iter
              (for nrow% :in nrows)
              (for object :in objects)
              (fill-matrix-elements-using% object elements 
                                           offset nrow nil)
              (incf offset nrow%))
            result)))))

(defun stack-horizontally (&rest objects)
  "Stack arguments horizontally, converting to a common type.  A
  vector is interpreted as a column matrix."
  (stack :h objects))

(defun stack-vertically (&rest objects)
  "Stack arguments horizontally, converting to a common type.  A
  vector is interpreted as a column matrix."
  (stack :v objects))

;;; identity

(defun eye (n &key (kind :dense) (initial-element 1)
            lla-type)
  "Return an identity matrix of the given KIND (can also be :diagonal)
and LLA-TYPE.  INITIAL-ELEMENT will be used for the diagonal."
  (if (eq kind :diagonal)
      (make-diagonal n lla-type initial-element)
      (bind (((:lla-matrix eye) (make-matrix n n lla-type :kind kind))
             (initial-element (coerce* initial-element lla-type)))
        (dotimes (i n)
          (setf (eye (eye-index i i)) initial-element))
        eye)))

;;; vector concatenation

(defun concat (&rest vectors)
  "Concatenate VECTORS in to a single vector, using the narrowest common element
type supported by the implementation.  Lists are treated as SIMPLE-VECTORS,
atoms that are not vectors wrapped in a vector."
  (let* ((vectors (mapcar (lambda (v)
                            (typecase v
                              (vector v)
                              (list (coerce v 'simple-vector))
                              (t (vector v))))
                          vectors))
         (lla-types (mapcar #'array-lla-type vectors))
         (lengths (mapcar #'length vectors))
         (common-type (reduce #'common-lla-type lla-types))
         (result (lla-vector (reduce #'+ lengths) common-type)))
    (iter
      (with offset := 0)
      (for v :in vectors)
      (for l :in lengths)
      (copy-elements v 0 result offset l)
      (incf offset l))
    result))

;;; grouping

(defgeneric group-by-index (object index &optional m)
  (:documentation "Group elements of object by index.  Return groups a
  simple-vector of similar objects.  M is the number of indices, which start at
  0.  Order is preserved within groups."))

(defmethod group-by-index ((matrix dense-matrix-like) index &optional
                              (m (1+ (emax index))))
  "Separate the rows of a matrix into submatrices, based on index."
  (let ((submatrices (make-array m :initial-element nil)))
    (dotimes (row-index (nrow matrix))
      (let* ((row (sub matrix row-index t)))
        (push row (aref submatrices (aref index row-index)))))
    (map 'vector (lambda (sm)
                   (stack :v (nreverse sm)))
         submatrices)))

(defmethod group-by-index ((vector vector) index &optional
                           (m (1+ (emax index))))
  (let ((subvectors (make-array m :initial-element nil))
        (result-type `(simple-array ,(array-element-type vector) (*))))
    (map nil (lambda (element i)
               (push element (aref subvectors i)))
         vector index)
    (map 'vector (lambda (subvector)
                   (coerce (nreverse subvector) result-type))
         subvectors)))
