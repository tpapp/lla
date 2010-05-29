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

(defun stack% (objects horizontal?)
  (bind (((:values nrows ncols)
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
                 (result (make-matrix lla-type nrow 
                                      (reduce #'+ ncols)))
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
                 (result (make-matrix lla-type nrow ncol))
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
  (stack% objects t))

(defun stack-vertically (&rest objects)
  "Stack arguments horizontally, converting to a common type.  A
  vector is interpreted as a column matrix."
  (stack% objects nil))

;;; identity

(defun eye (n &key (kind :dense) (initial-element 1)
            lla-type)
  "Return an identity matrix of the given KIND (can also be :diagonal)
and LLA-TYPE.  INITIAL-ELEMENT will be used for the diagonal."
  (if (eq kind :diagonal)
      (make-diagonal lla-type n initial-element)
      (bind (((:lla-matrix eye) (make-matrix lla-type n n :kind kind))
             (initial-element (coerce* initial-element lla-type)))
        (dotimes (i n)
          (setf (eye (eye-index i i)) initial-element))
        eye)))
