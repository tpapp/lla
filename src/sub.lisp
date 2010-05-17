;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla)

(defmethod sub ((nv numeric-vector) &rest ranges)
  (bind (((range) ranges)
         ((:slots-r/o lla-type elements) nv)
         (range (transform-range range (length elements))))
    (if (typep range 'fixnum)
        (aref elements range)
        (bind ((result-length (range-dimension range))
               ((:lla-vector result :elements result-elements) (make-nv lla-type result-length)))
          (etypecase range
            (sub-range (copy-elements result-length elements (car range) lla-type
                                      result-elements 0))
            (vector (iter
                      (for index :in-vector range)
                      (for result-index :from 0)
                      (setf (aref result-elements result-index)
                            (aref elements index)))))
          result))))

(defmethod sub ((matrix dense-matrix-like) &rest ranges)
  ;; we traverse the matrix as if it was transposed and row-major
  (declare (optimize debug))
  (bind (((row-range col-range) ranges)
         ((:slots-r/o lla-type elements nrow ncol) matrix))
    (with-range-indexing ((list col-range row-range) (vector ncol nrow)
                          matrix-inc matrix-index end? result-dimensions)
      (bind (((:values result-length result)
              (ecase (length result-dimensions)
                (0 (return-from sub (aref elements (matrix-index))))
                (1 (bind ((#(length) result-dimensions))
                     (values length (make-nv lla-type length))))
                (2 (bind ((#(result-ncol result-nrow) result-dimensions))
                     (values (* result-ncol result-nrow) (make-matrix lla-type result-nrow result-ncol
                                                                      :kind (matrix-kind matrix)))))))
             (result-elements (elements result)))
        (iter
          (for result-index :from 0 :below result-length)
          (setf (aref result-elements result-index) (aref elements (matrix-index)))
          (matrix-inc))
        result))))
