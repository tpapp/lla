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

(defun make-diagonal (length lla-type &optional initial-element)
  "Make a diagonal (usage of INITIAL-ELEMENT follows the semantics of
LLA-VECTOR)."
  (make-diagonal% (lla-vector length lla-type initial-element)))

(defmethod pack ((diagonal diagonal))
  (make-diagonal% (pack (elements diagonal))))
