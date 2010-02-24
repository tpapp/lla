;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla)

(defmethod metabang.bind.developer:bind-generate-bindings ((kind (eql :lla-vector))
                                                           variable-form value-form
                                                           body declarations remaining-bindings)
  (assert value-form () "No values provided!")
  (bind (((var) variable-form)
         (elements (gensym* 'elements- var)))
    (check-type var symbol)
    `((let* ((,var ,value-form)
              (,elements (elements ,var)))
         (macrolet ((,var (index)
                      `(aref ,',elements ,index)))
           ,(metabang-bind::bind-filter-declarations declarations variable-form)
           ,@(metabang-bind::bind-macro-helper remaining-bindings declarations body))))))

(defmethod metabang.bind.developer:bind-generate-bindings ((kind (eql :lla-matrix))
                                                           variable-form value-form
                                                           body declarations remaining-bindings)
  (assert value-form () "No values provided!")
  (bind (((var) variable-form)
         (elements (gensym* 'elements- var))
         (nrow (gensym* 'nrow- var)))
    (check-type var symbol)
    `((let* ((,var ,value-form)
             (,elements (elements ,var))
             (,nrow (nrow ,var)))
        (macrolet ((,var (row col)
                     `(aref ,',elements (cm-index2 ,',nrow ,row ,col))))
          ,(metabang-bind::bind-filter-declarations declarations variable-form)
           ,@(metabang-bind::bind-macro-helper remaining-bindings declarations body))))))
