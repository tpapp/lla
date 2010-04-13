;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla)

(setf (binding-form-docstring :lla-vector)
      (cons :lla-vector
            "(:lla-vector var &key elements accessor length) will bind
            the value form to VAR, its elements to ELEMENTS (if given,
            default is NIL), and bind a local function as the
            accessor (ACCESSOR INDEX) (default is VAR, unless ELEMENTS
            is given, when it is NIL)."))

(defmethod metabang.bind.developer:bind-generate-bindings 
    ((kind (eql :lla-vector))
     variable-form value-form
     body declarations remaining-bindings)
  (bind (((value) value-form)
         ((var &key elements
               (accessor (unless elements var))
               length) variable-form)
         (elements (aif elements
                        it
                        (gensym* 'elements- var))))
    (check-type var symbol)
    `((let* ((,var ,value)
             (,elements (elements ,var))
             ,@(when length
                 `((,length (length ,elements)))))
        (,@(when accessor
             `(flet ((,accessor (index)
                       (aref ,elements index))
                     ((setf ,accessor) (value index)
                       (setf (aref ,elements index) value)))
                (declare (ignorable (function ,accessor) (function (setf ,accessor))))))
           ,(metabang-bind::bind-filter-declarations declarations variable-form)
           ,@(metabang-bind::bind-macro-helper 
              remaining-bindings declarations body))))))

(defmethod metabang.bind.developer:bind-generate-bindings
    ((kind (eql :lla-matrix))
     variable-form value-form
     body declarations remaining-bindings)
  (bind (((value) value-form)
         ((var &key elements
               (accessor (if (not elements) var))
               (indexer (make-symbol* var '#:-index))
               length (nrow (gensym* '#:nrow- var))) variable-form)
         (elements (aif elements
                        it
                        (gensym* 'elements- var))))
    (check-type var symbol)
    `((let* ((,var ,value)
             (,elements (elements ,var))
             (,nrow (nrow ,var))
             ,@(if length
                   `((,length (length ,elements)))
                   nil))
        (flet (,@(when accessor
                   `((,accessor (index)
                                (aref ,elements index))
                     ((setf ,accessor) (value index)
                      (setf (aref ,elements index) value))))
               ,@(when indexer
                   `((,indexer (row col)
                               (cm-index2 ,nrow row col)))))
          ,@(when accessor
              `(declare (ignorable (function ,accessor) (function (setf ,accessor)))))
          ,@(when indexer
              `(declare (ignorable (function ,indexer))))
          ,(metabang-bind::bind-filter-declarations declarations variable-form)
          ,@(metabang-bind::bind-macro-helper 
             remaining-bindings declarations body))))))
