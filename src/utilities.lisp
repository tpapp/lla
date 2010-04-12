(in-package :lla)

(define-make-symbol% :lla)

(eval-when (:compile-toplevel :load-toplevel :execute)
  (pushnew :muffle-notes cl:*features*))

(defun relative-value% (value new-value relative-p)
  "Calculate new value, depending on RELATIVE-P."
  (if relative-p
      (+ value new-value)
      new-value))
