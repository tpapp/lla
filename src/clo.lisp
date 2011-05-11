;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla)

(defun remove-row-separator (elements &optional (separator :/))
  "Remove SEPARATOR from elements.  The two other dimensions are returned as the
  second and the third value.  If no separator was found, the third value is NIL and
  the second value is the length of elements."
  (let ((ncol (position separator elements))
        (length (length elements)))
    (if ncol
        (bind ((elements (remove separator elements))
               (length (1- length))
               ((:values nrow remainder) (floor length ncol)))
          (assert (plusp ncol) ()
                  "Zero number of columns.")
          (assert (not (find separator elements)) ()
                  "Multiple separators (~S) found." separator)
          (assert (zerop remainder) ()
                  "The number of elements (~A) is not divisible by NCOL (~A)."
                  length ncol)
          (values (coerce (remove separator elements) 'vector) nrow ncol))
        (values (coerce elements 'vector) (length elements)))))


(defun strip-keys (list key-lists &key (test #'eq) defaults)
  "Traverse list as long as keys are found in KEY-LISTs, then return the rest of the
list and a vector of keys, positions corresponding to KEY-LISTS.  Check that at most
one key occurs from each key list.  If no key occurs, DEFAULTS is used, should be a
sequence or a single element, the latter is recycled."
  ;; !! code is a bit ugly
  (let* ((n (length key-lists))
         (flags (make-array n :initial-element 0 :element-type 'bit))
         (keys (make-array n :initial-contents
                           (if (typep defaults '(and sequence (not null)))
                               defaults
                               (make-sequence 'list n :initial-element defaults)))))
    (loop
      while (and list
                 (not 
                   (iter
                     (for key-list :in key-lists)
                     (for position :from 0)
                     (let ((key (find (first list) key-list :test test)))
                       (when key
                         (assert (zerop (aref flags position)) ()
                                 "Multiple keys found: ~A and ~A both in ~A."
                                 (aref keys position) key key-list)
                         (setf (aref keys position) key
                               (aref flags position) 1))
                       (never key)))))
      do (setf list (cdr list)))
    (values list keys)))

;; (strip-keys '(1 2 3) nil)

;; (strip-keys '(d a 1 2 3)
;;             '((a b)
;;               (c d)))

;; (strip-keys '(a d b 1 2 3)
;;             '((a b)
;;               (c d)))

(defmacro clo (&rest arguments)
  "Create LLA Object.  Syntax:

  arguments ::= keyword* elements
  keyword ::= kind | lla-type
  elements ::= element* | element* :/ element*

  Possible KIND keywords
  are :diagonal, :dense, :lower(-triangular), :upper(-triangular) and :hermitian.  If
  none are given, the result is a numeric-vector, unless an element separator (:/) is
  given.  At most one KIND keyword can be given.

  At most one LLA-TYPE keyword can be given, if one is supplied, it
  elements will be coerced to that type.

  :/ in elements terminates a row for matrix types.  It can occur at
  most once.  If not given, all elements will be in a single row."
  (bind (((:values elements keys)
          (strip-keys arguments (load-time-value (list (lla-types)
                                                       '(:diagonal :dense
                                                         :lower :lower-triangular
                                                         :upper :upper-triangular
                                                         :hermitian)))
                      :defaults '(t :dense)))
         (#(element-type kind) keys)
         ((:values elements nrow ncol) (remove-row-separator elements))
         (dimensions (if ncol (list nrow ncol) nrow))
         (elements-var (gensym* '#:elements-var))
         represented-elements
         ((:flet save-element (index element))
          (push (list index element) represented-elements)))
    ;; collect elements that are set
    (if ncol
        (row-major-loop (dimensions index row col)
          (when (or (eql kind :dense) (represented-element? kind row col))
            (save-element index (aref elements index))))
        (loop for index :from 0
              for element :across elements
              do (save-element index element)))
    `(let ((,elements-var (make-array* ',dimensions ',element-type)))
       ;; fill elements
       ,@(loop
           for (index element) :in (nreverse represented-elements)
           collect `(setf (row-major-aref ,elements-var ,index)
                          (coerce* ,element ,element-type)))
       ;; create object !! will be used when introducing kinds
       ,(if (eq kind :dense)
            elements-var
            `(make-matrix ,kind ',dimensions :initial-contents ,elements-var)))))
