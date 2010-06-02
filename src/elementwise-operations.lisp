;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla)

(defgeneric emap-elements% (object)
  (:documentation "Return (values TYPE ELEMENTS DIMENSIONS) for
  objects which are analogous to vectors, or (values nil ELEMENTS nil)
  otherwise.  TYPE is a symbol, uniquely mapped to the type of OBJECT,
  ELEMENTS is a SIMPLE-ARRAY1, containing the elements of OBJECT.
  DIMENSION is a list of dimensions, or a single dimension (will be
  compared with EQUAL)")
  (:method (object)
    (values nil object nil))
  (:method ((vector vector))
    (values 'vector
            (as-simple-array1 vector)
            (array-dimension vector 0)))
  (:method ((array array))
    (values 'array
            (as-simple-array1 array)
            (array-dimensions array)))
  (:method ((diagonal diagonal))
    (bind (((:accessors-r/o elements) diagonal))
      (values 'diagonal
              elements
              (length elements))))
  (:method ((matrix dense-matrix-like))
    (set-restricted matrix)
    (bind (((:accessors-r/o elements nrow ncol) matrix))
      (values (type-of matrix)
              elements
              (list nrow ncol)))))

(defgeneric emap-dense-type% (function type scalars?)
  (:documentation "EMAP uses this function to recognize certain operations that
  preserve the matrix type.  It is called when TYPE is a subtype of
  DENSE-MATRIX-LIKE.  SCALARS? is nil only when all ARGs are matrices.")
  (:method (function type scalars?)
    ;; fallback method: always assume dense, unless have reason not to
    (declare (ignore function type scalars?))
    'dense-matrix)
  (:method ((function (eql #'+)) type scalars?)
    (if scalars?
        'dense-matrix
        type))
  (:method ((function (eql #'-)) type scalars?)
    (if scalars?
        'dense-matrix
        type))
  (:method ((function (eql #'*)) type scalars?)
    (declare (ignore scalars?))
    type)
  (:method ((function (eql #'/)) type scalars?)
    ;; the correctness of this depends on (/ ... 0) leading to an error
    (declare (ignore scalars?))
    type))

(defgeneric emap-result% (type elements dimensions)
  (:documentation "Create return value for EMAP.")
  (:method ((type (eql 'vector)) elements dimensions)
    elements)
  (:method ((type (eql 'array)) elements dimensions)
    (displace-array elements dimensions))

  (:method ((type (eql 'diagonal)) elements dimensions)
    (make-diagonal% elements))
  (:method ((type (eql nil)) elements dimensions)
    elements))

(defmacro define-dense-emap-result% (kind)
  `(defmethod emap-result% ((type (eql ',(matrix-type kind)))
                            elements dimensions)
     (bind (((nrow ncol) dimensions))
       (make-matrix% nrow ncol elements :kind ,kind))))

(define-dense-emap-result% :dense)
(define-dense-emap-result% :hermitian)
(define-dense-emap-result% :lower)
(define-dense-emap-result% :upper)

(defun emap-vectors% (function args length)
  "Map ARGS into a simple vector elementwise, using FUNCTION.  The
element type of the result the narrowest common (extended) LLA type
that can accommodate all results.  Some of the args can be atoms,
in which case they are treated as a vector filled with that atom."
  (declare (optimize speed))
  (if length
      (let (lla-type result)
        (declare (fixnum length))
        (iter
          (for index :from 0 :below length)
          (declare (iterate:declare-variables)
                   (fixnum index))
          (let ((element 
                 (apply function 
                        (mapcar (lambda (arg) 
                                  (if (simple-array? arg)
                                      (row-major-aref arg index)
                                      arg))
                                args))))
            ;; when first time, create array
            (when (zerop index)
              (setf lla-type (atom-representable-lla-type element)
                    result (lla-vector lla-type length)))
            (locally
                (declare (type (simple-array * (*)) result))
                (handler-case (setf (row-major-aref result index)
                                    (coerce* element lla-type))
                  (type-error ()
                    (let ((old-result result))
                      (setf lla-type 
                            (common-lla-type 
                             lla-type 
                             (atom-representable-lla-type element))
                            result (lla-vector lla-type length))
                      (dotimes (i index)
                        (setf (aref result i) 
                              (coerce* (aref old-result i) lla-type)))
                      (setf (aref result index) element)))))))
        result)
      (apply function args)))

(defun emap-common-type% (common-type type)
  "Find the accumulated common type, if possible, otherwise signal an error."
  (cond
    ((null type) common-type)           ; scalar, ignored
    ((null common-type) type)           ; all scalars so far
    ((and (subtypep common-type 'dense-matrix-like) 
          (subtypep type 'dense-matrix-like))
     (if (eq common-type type)
         type
         'dense-matrix))
    ((eq common-type type) type)
    (t "Type ~A is not compatible with ~A." type common-type)))

(defun emap (function &rest args)
  "Map ARGS into a similar object elementwise, using FUNCTION.  The element type
of the result the narrowest common (extended) LLA type that can accommodate all
results.  Some of the args can be atoms, in which case they are treated as an
object filled with that atom.  All ARGs need to have the same type, see
EMAP-ELEMENTS%.  FORCE-DENSE? "
  (let* (type
         dimensions
         scalars?
         (elements
          (iter
            (for arg :in args)
            (bind (((:values arg-type arg-elements arg-dimensions)
                    (emap-elements% arg)))
              (if arg-type
                  (if type
                      (assert (equal arg-dimensions dimensions) ()
                              "Dimensions ~A are not compatible ~
                                 with ~A" arg-dimensions dimensions)
                      (setf dimensions arg-dimensions))
                  (setf scalars? t))
              (setf type (emap-common-type% type arg-type))

              (collecting arg-elements))))
         (result (emap-vectors% function elements
                                (cond 
                                  ((null type)
                                   nil)
                                  ((listp dimensions)
                                   (reduce #'* dimensions))
                                  (t dimensions)))))
    (when (subtypep type 'dense-matrix-like)
      (setf type (emap-dense-type% function type scalars?)))
    (emap-result% type result dimensions)))

(defmacro define-elementwise-operation (name arglist documentation
                                        op-form)
  `(defgeneric ,name ,arglist
     (:documentation ,documentation)
     (:method ,arglist
       ,op-form)))

;;; shorthand for specific elementwise operations


(define-elementwise-operation e+ (object &rest objects)
  "Elementwise +."
  (apply #'emap #'+ object objects))

(define-elementwise-operation e- (object &rest objects)
  "Elementwise -."
  (apply #'emap #'- object objects))

(define-elementwise-operation e* (object &rest objects)
  "Elementwise *."
  (apply #'emap #'* object objects))

(define-elementwise-operation e/ (object &rest objects)
  "Elementwise /."
  (apply #'emap #'/ object objects))

(define-elementwise-operation eexpt (base power)
  "Elementwise EXPT."
  (emap #'expt base power))

(define-elementwise-operation eexp (arg)
  "Elementwise EXP."
  (emap #'exp arg))

(define-elementwise-operation elog (arg)
  "Elementwise LOG."
  (emap #'log arg))

;;; !!! write optimized 2-argument versions

;;; ?? Theoretically, a lower triangular matrix multiplied by an upper
;;; triangular one would be diagonal, but here constrain ourselves to returning
;;; a dense-matrix-like object. ??? Maybe we could return a diagonal here, I
;;; haven't fixed this in x
