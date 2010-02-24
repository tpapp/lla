;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

;;; Conversion and creation methods for the XARRAY interface.

(in-package #:lla)

(defun parse-nv-options (options &optional (default :double))
  "Parse XCREATE options for the NUMERIC-VECTOR class."
  (bind (((&key (lla-type default)) options))
    lla-type))

(defun xsimilar% (rank object)
  ;; not exported
  (let ((lla-type (lla-type object))
        (rank (if (eq t rank) (xrank object) rank)))
    (ecase rank
      (1 (list (nv-class lla-type)))
      (2 `(dense-matrix :lla-type ,lla-type)))))

(defmethod xsimilar (rank (nv numeric-vector))
  (xsimilar% rank nv))

(defmethod xsimilar (rank (m dense-matrix-like))
  (xsimilar% rank m))

(defmethod xcreate ((class (eql 'numeric-vector)) dimensions &optional
                    options)
  (bind (((n) dimensions))
    (make-nv n (parse-nv-options options))))

(defmethod xcreate ((class (eql 'dense-matrix)) dimensions &optional options)
  (unless (= (length dimensions) 2)
    (error "Exactly 2 dimensions are needed for a matrix."))
  (bind (((nrow ncol) dimensions)
         ((&key (lla-type :double)) options))
    (make-matrix  nrow ncol lla-type)))

(defmethod as* ((class (eql 'numeric-vector)) 
                (object dense-matrix-like) copy-p options)
  (set-restricted object)
  (copy-nv object :destination-type (parse-nv-options options (lla-type object))
           :copy-p copy-p))

;;; diagonal

(defmethod xcreate ((class (eql 'diagonal)) dimensions &optional
                    options)
  (bind (((n) dimensions))
    (make-diagonal n (parse-nv-options options))))

(defmethod xsimilar ((rank (eql 1)) (object diagonal))
  `(diagonal :lla-type ,(lla-type object)))

(defmethod as* ((class (eql 'numeric-vector)) (object diagonal) copy-p options)
  (copy-nv object :destination-type (parse-nv-options options (lla-type object))
           :copy-p copy-p))

;;; matrix -- !!! needs to be written
