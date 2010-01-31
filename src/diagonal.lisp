(in-package :lla)

;;; Diagonals are not matrices: they just *represent* a diagonal of a
;;; matrix, or a diagonal matrix.  They are a subclass of
;;; NUMERIC-VECTOR.

(define-abstract-class diagonal (numeric-vector)
  ()
  (:documentation "Elements representing a diagonal of a matrix or a
  diagonal matrix."))

(define-lla-class diagonal)

(defun diagonal-class (lla-type)
  "Return the class name corresponding to LLA-TYPE."
  (append-lla-type diagonal lla-type))

(defun nv->diagonal% (nv)
  "Convert numeric vector to diagonal.  Not exported, internal use
only.  Destructive, changes class of NV."
  (change-class nv (diagonal-class (lla-type nv))))

(defun make-diagonal* (lla-type elements)
  "Like MAKE-NV*, but for diagonals."
  (make-instance (diagonal-class lla-type) :elements elements))

(defmethod make-load-form ((diagonal diagonal) &optional environment)
  (declare (ignore environment))
  `(make-diagonal* ,(lla-type diagonal) ,(make-elements-load-form diagonal)))

(defun nv->diagonal (nv &key (copy-p nil))
  "Convert numeric vector to diagonal."
  (make-diagonal* (lla-type nv) (if copy-p
                                    (copy-elements nv)
                                    (elements nv))))
(defun diagonal->nv (diagonal &key (copy-p nil))
  "Convert numeric vector to diagonal."
  (make-nv* (lla-type diagonal) (if copy-p
                              (copy-elements diagonal)
                              (elements diagonal))))

(defun make-diagonal (length lla-type &optional (initial-element 0))
  "Make a diagonal."
  (nv->diagonal% (make-nv length lla-type initial-element)))

(defun create-diagonal (initial-contents &optional lla-type)
  "Create a diagonal from initial-contents."
  (nv->diagonal% (create-nv initial-contents lla-type)))
                
(defun diagonal->matrix (diagonal &optional (kind :dense))
  "Create a matrix of given KIND from DIAGONAL."
  (let* ((diagonal-elements (elements diagonal))
         (n (length diagonal-elements))
         (matrix (make-matrix n n (lla-type diagonal) :kind kind))
         (matrix-elements (elements matrix)))
    (dotimes (i n)
      (setf (aref matrix-elements (cm-index2 n i i))
            (aref diagonal-elements i)))
    matrix))

(defun matrix->diagonal (matrix)
  "Extract the diagonal from a matrix."
  (bind (((:slots-read-only nrow ncol (matrix-elements elements)) matrix)
         (n (min nrow ncol))
         (lla-type (lla-type matrix))
         (diagonal-elements (make-nv-elements n lla-type)))
    (dotimes (i n)
      (setf (aref diagonal-elements i)
            (aref matrix-elements (cm-index2 n i i))))
    (make-diagonal* lla-type diagonal-elements)))

(defmethod xcreate ((class (eql 'diagonal)) dimensions &optional
                    options)
  (bind (((n) dimensions))
    (make-diagonal n (extract-only-lla-type options))))

(defmethod xsimilar ((rank (eql 1)) (object diagonal))
  `(diagonal :lla-type ,(lla-type object)))
