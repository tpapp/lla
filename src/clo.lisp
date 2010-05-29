;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla)

(defun remove-row-separators% (elements &optional (separator :/))
  "Remove SEPARATORs from elements.  The position of the first
separator (or if none, the length of elements) is returned as the
second value."
  (let ((ncol (aif (position separator elements)
                   it
                   (length elements))))
    (assert (plusp ncol) ()
            "Zero number of columns.")
    (assert (<= (count separator elements) 1) ()
            "Multiple separators (~S) found." separator)
    (values (remove separator elements) ncol)))

(defun rearrange-elements% (elements ncol kind)
  "Rearrange ELEMENTS into a matrix.  Check that NCOL divides length
of ELEMENTS, return (values ELEMENTS NROW).  ELEMENTS is a list."
  (bind (((:accessors-r/o length) elements)
         ((:values nrow remainder) (floor length ncol))
         (row 0)
         (col 0)
         (result (make-array length)))
    (unless (zerop remainder)
      (error "The number of elements (~A) is not divisible by NCOL (~A)."
             length ncol))
    (iter
      (for element :in elements)
      (setf (aref result (cm-index2 nrow row col)) 
            (if (matrix-mask kind row col)
                element
                0))
      (incf col)
      (when (= col ncol)
        (incf row)
        (setf col 0)))
    (values (coerce result 'list) nrow)))

(defmacro clo (&rest arguments)
  "Create LLA Object.  Syntax:

  arguments ::= keyword* elements
  keyword ::= kind | lla-type | no-coerce
  elements ::= element* | element* :/ element*

  Possible KIND keywords are :diagonal, :dense, :lower-triangular
  and :upper-triangular.  If none are given, the result is a
  numeric-vector, unless an element separator (:/) is given.  At most
  one KIND keyword can be given.

  At most one LLA-TYPE keyword can be given, if one is supplied, it
  elements will be coerced to that type.

  :/ in elements terminates a row for matrix types.  It can occur at
  most once.  If not given, all elements will be in a single row."
  (macrolet ((set-if-first (place arg &key (name (symbol-name place))
                                  single-value?)
               (check-type place symbol)
               `(if ,place
                    ,(if single-value?
                         `(error "Multiple ~S declarations." ,arg)
                         `(error "Multiple ~A declarations (~S,~S)" 
                                 ,name ,place ,arg))
                    (setf ,place ,arg))))
    (bind (kind lla-type
           (elements
            (iter
              (for on-arg :on arguments)
              (for arg := (car on-arg))
              (cond
                ((not (keywordp arg))
                 (return on-arg))
                ((eq :diagonal arg)
                 (set-if-first kind arg))
                ((valid-matrix-kind? arg)
                 (set-if-first kind arg))
                ((typep arg 'lla-type)
                 (set-if-first lla-type arg))
                (t (error "Unrecognized keyword: ~A" arg)))))
           (elements-var (gensym* '#:elements-var))
           ((:flet elements-form (elements creation-form))
            (assert (notany #'keywordp elements) ()
                    "Elements ~A contain keywords." elements)
             `(let* ((,elements-var 
                      ,(if lla-type
                           `(map '(simple-array ,(lla->lisp-type lla-type) (*))
                                 (lambda (x)
                                   (coerce* x ,lla-type))
                                 (vector ,@elements))
                           `(vector ,@elements))))
                ,creation-form)))
      ;; no kind specification, but separator found => dense matrix
      (when (and (not kind) (find :/ elements))
        (setf kind :dense))
      ;; return form
      (case kind
        ((nil)
           (elements-form elements elements-var))
        (:diagonal
           (elements-form elements `(make-diagonal% ,elements-var)))
        (otherwise
           (bind (((:values elements ncol)
                   (remove-row-separators% elements))
                  ((:values elements nrow)
                   (rearrange-elements% elements ncol kind)))
             (elements-form elements
                            `(make-matrix% ,nrow ,ncol
                                           ,elements-var :kind ,kind))))))))
