;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla)



(setf (binding-form-docstring :lla-vector)
      (cons :lla-vector
            "(:lla-vector var &key accessor length) will bind the
            value form to VAR and bind a local function as the
            accessor (ACCESSOR INDEX) (default is VAR, unless ELEMENTS
            is given, when it is NIL)."))

(defmethod metabang.bind.developer:bind-generate-bindings 
  ((kind (eql :lla-vector))
   variable-form value-form
   body declarations remaining-bindings)
  (bind (((value) value-form)
         ((var &key
               (accessor var)
               length
               declarations?) variable-form))
    (check-type var symbol)
    `((let* ((,var ,value)
             ,@(when length
                 `((,length (length ,var)))))
        (flet (,@(when accessor
                   `((,accessor (index)
                                (aref ,var index))
                     ((setf ,accessor) (value index)
                      (setf (aref ,var index) value)))))
          ,@(when accessor
              `((declare (ignorable (function ,accessor) (function (setf ,accessor))))))
          ,(metabang-bind::bind-filter-declarations declarations variable-form)
          ,@(maybe-wrap-list
             (when declarations? `(with-vector-type-declarations (,var)))
             (metabang-bind::bind-macro-helper 
              remaining-bindings declarations body)))))))

(defmethod metabang.bind.developer:bind-generate-bindings
  ((kind (eql :lla-matrix))
   variable-form value-form
   body declarations remaining-bindings)
  (bind (((value) value-form)
         ((var &key
               (elements (gensym* 'elements- var) elements?)
               (accessor (unless elements? var))
               (indexer (make-symbol* var '#:-index))
               length (nrow (gensym* '#:nrow- var))
               ncol
               declarations?) variable-form)
         (elements (aif elements
                        it
                        (gensym* 'elements- var))))
    (check-type var symbol)
    `((let* ((,var ,value)
             (,elements (elements ,var))
             (,nrow (nrow ,var))
             ,@(if length
                   `((,length (length ,elements)))
                   nil)
             ,@(if ncol
                   `((,ncol (ncol ,var)))))
        (flet (,@(when accessor
                   `((,accessor (index)
                                (aref ,elements index))
                     ((setf ,accessor) (value index)
                      (setf (aref ,elements index) value))))
               ,@(when indexer
                   `((,indexer (row col)
                               (cm-index2 ,nrow row col)))))
          ,@(when accessor
              `((declare (ignorable (function ,accessor) (function (setf ,accessor))))))
          ,@(when indexer
              `((declare (ignorable (function ,indexer)))))
          ,@(when ncol
              `((declare (ignorable ncol))))
          ,(metabang-bind::bind-filter-declarations declarations variable-form)
          ,@(maybe-wrap-list
             (when declarations? `(with-vector-type-declarations (,elements)))
             (metabang-bind::bind-macro-helper 
              remaining-bindings declarations body)))))))
