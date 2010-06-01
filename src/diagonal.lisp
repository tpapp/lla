(in-package :lla)

;;; Diagonals are not matrices: they just *represent* a diagonal of a
;;; matrix, or a diagonal matrix.  They are a subclass of
;;; NUMERIC-VECTOR-LIKE.

(defstruct (diagonal
             (:conc-name)
             (:constructor make-diagonal% (elements))
             (:copier nil)
             (:include elements%)
             (:print-object))
  "Object representing a diagonal of a matrix or a diagonal matrix.")

(defun make-diagonal (lla-type length &optional (initial-element 0))
  "Make a diagonal."
  (make-diagonal% (lla-vector lla-type length initial-element)))

(defmethod pack ((diagonal diagonal))
  (make-diagonal% (pack (elements diagonal))))
