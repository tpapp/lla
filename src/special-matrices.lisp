;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package :lla)

;;;; mean accumulator

;; (defstruct (wrapped-matrix-mean-accumulator
;;              (:constructor wrapped-matrix-mean-accumulator% (mean kind))
;;              (:include array-mean-accumulator))
;;     "Accumulator for wrapped matrices.  Save type as extra information."
;;   (kind nil))

;; (defun wrapped-matrix-mean-accumulator (wrapped-matrix)
;;   "Create a conforming wrapped-matrix-mean-accumulator."
;;   (wrapped-matrix-mean-accumulator%
;;    (make-array (array-dimensions (wrapped-matrix-elements wrapped-matrix))
;;                :initial-element 0d0)
;;    (matrix-kind wrapped-matrix)))

;; (define-conforming-accumulator (mean (matrix wrapped-matrix))
;;   (wrapped-matrix-mean-accumulator matrix))

;; (defmethod add :after ((accumulator wrapped-matrix-mean-accumulator) object)
;;   (let+ (((&structure wrapped-matrix-mean-accumulator- kind) accumulator))
;;     (when (and kind (not (eq kind (matrix-kind object))))
;;       (setf kind nil))))

;; (defmethod mean ((accumulator wrapped-matrix-mean-accumulator))
;;   (let+ (((&structure wrapped-matrix-mean-accumulator- mean kind) accumulator))
;;     (if kind
;;         (convert-matrix kind mean)
;;         mean)))

;;;; matrix creation and element access

(defgeneric matrix-kind (matrix)
  (:method ((matrix array))
    (check-type matrix aops:matrix)
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

(defgeneric convert-matrix (kind object)
  (:documentation "Convert object to a matrix of the given KIND.  May share
  structure unless the second value is T.")
  (:method ((kind (eql 'dense)) object)
    (aops:as-array object)))

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



;;;; special matrix types



;;;; Triangular matrices





;;;; Convenience functions for creating matrices

(defmacro dense (element-type &body rows)
  "Macro for creating a (dense) matrix (ie a rank 2 array).  ROWS should be a
list of lists (or atoms, which are treated as lists), elements are evaluated."
  (let+ ((rows (map 'vector #'ensure-list rows))
         (nrow (length rows))
         (ncol (length (aref rows 0)))
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

(defun diag (element-type &rest elements)
  "Return a DIAGONAL with elements coerced to ELEMENT-TYPE."
  (make-diagonal (apply #'vec element-type elements)))
