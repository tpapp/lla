;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla)

;;; extension of some generic functions in CL-numlib, also overriding
;;; for vectors

(defmethod size ((matrix dense-matrix-like))
  (* (nrow matrix) (ncol matrix)))

(defmethod size ((diagonal diagonal))
  (length (elements diagonal)))

(defun sum-elements% (vector start end)
  (declare (optimize speed (safety 0))
           (fixnum start end))
  (with-vector-type-expansion (vector)
    (lambda (lla-type)
      `(bind ((sum (zero* ,lla-type)))
         ,@(when lla-type
             `((declare (type ,(lla->lisp-type lla-type) sum))))
         (iter
           (for (the fixnum index) :from start :below end)
           (declare (iterate:declare-variables))
           (incf sum (aref vector index)))
         (muffle-optimization-notes
           sum)))))

(defmethod sum ((vector vector))
  (sum-elements% vector 0 (length vector)))

(defmethod sum ((matrix dense-matrix-like))
  (set-restricted matrix)
  (sum (elements matrix)))

(defmethod sum ((diagonal diagonal))
  (sum (elements diagonal)))

(eval-when (:compile-toplevel :load-toplevel :execute)
  )

(defun sse-elements% (vector start end &optional (mean 0))
  (declare (optimize speed (safety 0))
           (fixnum start end))
  (with-vector-type-expansion (vector)
    (lambda (lla-type)
      (let ((real-lla-type (real-lla-type lla-type)))
        `(bind ((sum (zero* ,real-lla-type))
                (mean (coerce* mean ,lla-type)))
           ,@(when lla-type
               `((declare (type ,(lla->lisp-type lla-type) mean))
                 (declare (type ,(lla->lisp-type real-lla-type) sum))))
           (iter
             (for (the fixnum index) :from start :below end)
             (declare (iterate:declare-variables))
             ,@(maybe-wrap
                (when (eq lla-type :integer)
                  '(muffle-optimization-notes))
                `(incf sum (let ((x (- (aref vector index) mean)))
                             ,(if (lla-complex? lla-type)
                                  ;; here we need the square *modulus*
                                  `(+ (expt (realpart x) 2)
                                      (expt (imagpart x) 2))
                                  `(expt x 2))))))
           (muffle-optimization-notes
             sum))))))

(defmethod sse ((vector vector) &optional (mean (mean vector)))
  (sse-elements% vector 0 (length vector) mean))

(defmethod sse ((matrix dense-matrix-like) &optional (mean (mean matrix)))
  (set-restricted matrix)
  (sse (elements matrix) mean))

(defmethod sse ((diagonal diagonal) &optional (mean (mean diagonal)))
  (sse (elements diagonal) mean))

