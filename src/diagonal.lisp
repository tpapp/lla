(in-package :lla)

;;; Diagonals are not matrices: they just *represent* a diagonal of a
;;; matrix, or a diagonal matrix.  They are a subclass of
;;; NUMERIC-VECTOR-LIKE.

(defclass diagonal (numeric-vector-like)
  ()
  (:documentation "Elements representing a diagonal of a matrix or a
  diagonal matrix."))

(defun make-diagonal* (lla-type elements)
  "Like MAKE-NV*, but for diagonals."
  (make-instance 'diagonal :lla-type lla-type :elements elements))

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

(defun make-diagonal (lla-type length &optional (initial-element 0))
  "Make a diagonal."
  (change-class (make-nv lla-type length initial-element)
                'diagonal))

(defun create-diagonal (initial-contents &optional lla-type)
  "Create a diagonal from initial-contents."
  (change-class (create-nv initial-contents lla-type)
                'diagonal))
                
(defun diagonal->matrix (diagonal &optional (kind :dense))
  "Create a matrix of given KIND from DIAGONAL."
  (let* ((diagonal-elements (elements diagonal))
         (n (length diagonal-elements))
         (matrix (make-matrix (lla-type diagonal) n n :kind kind))
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
         (diagonal-elements (make-nv-elements lla-type n)))
    (dotimes (i n)
      (setf (aref diagonal-elements i)
            (aref matrix-elements (cm-index2 n i i))))
    (make-diagonal* lla-type diagonal-elements)))
