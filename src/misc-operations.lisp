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

(defun map1-elements% (object function lla-type)
  "Map the elements using function, with the resulting LLA type."
  (map (nv-array-type lla-type) function (elements object)))

(defgeneric map1 (object function &optional lla-type)
  (:documentation "Map the object elementwise, returning an object of
  the same kind.  The lla-type of the result is the same as of the
  object by default.")
  (:method ((nv numeric-vector) function &optional (lla-type (lla-type nv)))
    (make-nv* lla-type (map1-elements% nv function lla-type)))
  (:method ((diagonal diagonal) function &optional (lla-type (lla-type diagonal)))
    (make-diagonal* lla-type (map1-elements% diagonal function lla-type)))
  (:method ((matrix dense-matrix-like) function &optional (lla-type (lla-type matrix)))
    (make-matrix* lla-type (nrow matrix) (ncol matrix)
                  (map1-elements% matrix function lla-type) :kind (matrix-kind matrix))))
