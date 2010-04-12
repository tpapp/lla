;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla)

(defmethod metabang.bind.developer:bind-generate-bindings 
    ((kind (eql :lla-vector))
     variable-form value-form
     body declarations remaining-bindings)
  (bind (((value) value-form)
         ((var &key elements
               (accessor (if (not elements) var))
               length) variable-form)
         (elements (aif elements
                        it
                        (gensym* 'elements- var))))
    (check-type var symbol)
    `((let* ((,var ,value)
             (,elements (elements ,var))
             ,@(if length
                   `((,length (length ,elements)))
                   nil))
        (,@(if accessor
               `(flet ((,accessor (index)
                         (aref ,elements index))
                       ((setf ,accessor) (value index)
                         (setf (aref ,elements index) value))))
               nil)
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
        (flet (,@(if accessor
                     `((,accessor (index)
                                  (aref ,elements index))
                       ((setf ,accessor) (value index)
                        (setf (aref ,elements index) value)))
                     nil)
               ,@(if indexer
                     `((,indexer (row col)
                                 (cm-index2 ,nrow row col)))))
           ,(metabang-bind::bind-filter-declarations declarations variable-form)
           ,@(metabang-bind::bind-macro-helper 
              remaining-bindings declarations body))))))
