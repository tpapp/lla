;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla)

(defmethod size ((nv numeric-vector-like))
  (length (elements nv)))

(defmethod sum ((nv numeric-vector-like))
  (declare (optimize speed (safety 0)))
  (bind (((:slots-r/o elements lla-type) nv))
    (expand-for-lla-types (lla-type :prologue (ecase lla-type))
      `(,lla-type
        (locally 
            (declare (type ,(nv-array-type lla-type) elements))
          (bind (((:accessors-r/o length) elements)
                 (sum (zero* ,lla-type)))
            (declare (,(lla-type->lisp-type lla-type) sum))
            (dotimes (i length)
              (incf sum (aref elements i)))
            sum))))))

(defmethod sse ((nv numeric-vector-like) &optional (mean (mean nv)))
  (declare (optimize speed (safety 0)))
  (bind (((:slots-r/o elements lla-type) nv))
    (expand-for-lla-types (lla-type :prologue (ecase lla-type))
      `(,lla-type
        (locally 
            (declare (type ,(nv-array-type lla-type) elements))
          (bind (((:accessors-r/o length) elements)
                 (sum (zero* ,lla-type))
                 (mean (coerce* mean ,lla-type)))
            (declare (,(lla-type->lisp-type lla-type) sum mean))
            (dotimes (i length)
              (incf sum (expt (- (aref elements i) mean) 2)))
            sum))))))
