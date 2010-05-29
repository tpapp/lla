;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla)

(defun emap-preprocess (args type dimensions-reader)
  "Preprocess args, comparing dimension for args of TYPE, return
extracted elements and the dimensions as the second value (which is
NIL, if no ARGS are of TYPE)."
  (let (dimensions)
    ;; determine and check length
    (values 
      (mapcar
       (lambda (arg)
         (if (typep arg type)
             (let ((arg-dimensions (funcall dimensions-reader arg)))
               (if dimensions
                   (assert (equal arg-dimensions dimensions)
                           () "Dimension mismatch for ~A" arg)
                   (setf dimensions arg-dimensions))
               (elements arg))))
       args)
      dimensions)))

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
      (values 'dense-matrix
              elements
              (list nrow ncol)))))

(defgeneric emap-result% (type elements dimensions)
  (:documentation "Create return value for EMAP.")
  (:method ((type (eql 'vector)) elements dimensions)
    elements)
  (:method ((type (eql 'array)) elements dimensions)
    (displace-array elements dimensions))
  (:method ((type (eql 'dense-matrix)) elements dimensions)
    (bind (((nrow ncol) dimensions))
      (make-matrix% nrow ncol elements)))
  (:method ((type (eql 'diagonal)) elements dimensions)
    (make-diagonal% elements))
  (:method ((type (eql nil)) elements dimensions)
    elements))

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
              (setf lla-type (representable-lla-type element)
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
                             (representable-lla-type element))
                            result (lla-vector lla-type length))
                      (dotimes (i index)
                        (setf (aref result i) 
                              (coerce* (aref old-result i) lla-type)))
                      (setf (aref result index) element)))))))
        result)
      (apply function args)))

(defun emap (function &rest args)
  "Map ARGS into a similar object elementwise, using FUNCTION.  The
element type of the result the narrowest common (extended) LLA type
that can accommodate all results.  Some of the args can be atoms, in
which case they are treated as an object filled with that atom.  All
ARGs need to have the same type, see EMAP-ELEMENTS%."
  (let* (type
         dimensions
         (elements
          (iter
            (for arg :in args)
            (bind (((:values arg-type arg-elements arg-dimensions)
                    (emap-elements% arg)))
              (when arg-type
                (if type
                    (progn
                      (assert (eq arg-type type) ()
                              "Type ~A is not compatible with ~A."
                              arg-type type)
                      (assert (equal arg-dimensions dimensions) ()
                              "Dimensions ~A are not compatible ~
                                 with ~A" arg-dimensions dimensions))
                    (setf type arg-type
                          dimensions arg-dimensions)))
              (collecting arg-elements))))
         (result (emap-vectors% function elements
                                (cond 
                                  ((null type)
                                   nil)
                                  ((listp dimensions)
                                   (reduce #'* dimensions))
                                  (t dimensions)))))
    (emap-result% type result dimensions)))

(defmacro define-elementwise-operation (name arglist documentation
                                        op-form)
  `(defgeneric ,name ,arglist
     (:documentation ,documentation)
     (:method ,arglist
       ,op-form)))

;;; shorthand for specific elementwise operations

;;; !!! TODO: have the result inherit matrix class whenever possible
;;; - Matrix-scalar operations:
;;; * and / ALWAYS preserve matrix kind (dense, lower/upper, hermitian)
;;; + and - NEVER preserve matrix kind, result is dense
;;;
;;; - Elementwise matrix operations:
;;; if matrices are of the same kind, so is the result
;;; if they are different kind, the result is dense, except for
;;; multiplication by lower/upper triangular matrices
;;;
;;; - Theoretically, a lower triangular matrix multiplied by an upper
;;; triangular one would be diagonal, but here constrain ourselves
;;; to returning a dense-matrix-like object. ??? Maybe we could
;;; return a diagonal here, I haven't fixed this in x

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



