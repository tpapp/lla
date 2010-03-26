(in-package :lla)

;;;; Printing and formatting

;;;; General variables and utility functions 

(defun print-length-truncate (dimension)
  "Return values (min dimension *print-length*) and whether the
constraint is binding."
  (if (or (not *print-length*) (<= dimension *print-length*))
      (values dimension nil)
      (values *print-length* t)))

(defvar *print-lla-precision* 5
  "number of digits after the decimal point when printing numeric matrices")

(defun standard-numeric-formatter (x)
  "Standard formatter for matrix printing.  Respects
*print-lla-precision*, and formats complex numbers as a+bi, eg
0.0+1.0i."
  ;; ?? do we want a complex numbers to be aligned on the +, like R? I
  ;; am not sure I like that very much, and for a lot of data, I would
  ;; visualize it graphically anyhow (I hate tables of 7+ numbers in
  ;; general).  -- Tamas, 2009-sep-13
  (typecase x
    (integer (format nil "~d" x))
    (real (format nil "~,vf" *print-lla-precision* x))
    (complex (format nil "~,vf+~,vfi"
		     *print-lla-precision* (realpart x)
		     *print-lla-precision* (imagpart x)))
    (t (format nil "~a" x))))


;;;; numeric vectors

(defun print-nv-elements (elements length stream)
  "Print elements of vector, automatically truncating to *PRINT-LENGTH*."
  (bind (((:values truncated-length truncated?) (print-length-truncate length)))
    (dotimes (i truncated-length)
      (princ (standard-numeric-formatter (aref elements i)) stream)
      (when (< (1+ i) truncated-length)
        (princ #\space stream)))
    (when truncated?
      (princ " ..." stream))))

(defmethod print-object ((obj numeric-vector) stream)
  (print-unreadable-object (obj stream :type t)
    (with-slots (elements) obj
      (let* ((length (length elements)))
	(format stream "~s with ~a elements: " (lla-type obj) length)
        (print-nv-elements elements length stream)))))

;;;; matrices

;;; None of the code below assumes anything about classes, it just
;;; uses the xref interface to access elements.

(defvar *print-matrix-aligned* t
  "If non-nil, characters will be aligned.")

(defvar *print-matrix-paddig* 1
  "Number of spaces between columns.")

(defun upper-triangular-mask (row col)
  (if (<= row col)
      nil
      "."))

(defun lower-triangular-mask (row col)
  (if (>= row col)
      nil
      "."))

(defun hermitian-mask (row col)
  ;; print * instead of . in the lower triangle
  (if (<= row col)
      nil
      "*"))

(defun print-matrix (matrix stream &key 
		     (formatter #'standard-numeric-formatter)
                     (accessor #'xref)
                     (mask (constantly nil)))
  "Format and print the elements of matrix (which is rank 2 xrefable
object, nothing further is needed) to stream, using
*print-matrix-padding* spaces between columns.  If
*print-matrix-aligned*, columns will be right-aligned.  Prints at most
*print-length* rows and columns, indicating more with a ...  Mask is a
function called with row and col indexes, and should return nil for
elements printed.  If mask returns a non-nil value, that will be
printed instead (should be a string)."
  ;; ?? maybe column & row labels, not a high priority at the moment
  (bind (((:values nrow row-trunc-p) (print-length-truncate (nrow matrix)))
	 ((:values ncol col-trunc-p) (print-length-truncate (ncol matrix)))
	 (formatted-elements (make-array (list nrow ncol)))
	 (column-widths (make-array ncol :element-type 'fixnum :initial-element 0))
	 (padding (make-array *print-matrix-paddig*
                              :element-type 'character :initial-element #\space))
	 (aligned-p *print-matrix-aligned*))
    ;; first pass - format elements, measure width
    (dotimes (col ncol)
      (dotimes (row nrow)
	(let* ((mask (funcall mask row col))
               (formatted-element (if mask
                                      mask
                                      (funcall formatter 
                                               (funcall accessor matrix row col))))
	       (width (length formatted-element)))
	  (setf (aref column-widths col) (max (aref column-widths col) width)
		(aref formatted-elements row col) formatted-element))))
    ;; second pass - print
    (dotimes (row nrow)
      (when (plusp row)
	(fresh-line stream))
      (dotimes (col ncol)
	(when (plusp col)
	  (princ padding stream))
	(let ((elt (aref formatted-elements row col)))
	  (if aligned-p
	      (format stream "~V@A" (aref column-widths col) elt)
	      (princ elt stream))))
      (when col-trunc-p
	(princ padding stream)
	(princ "..." stream)))
    (when row-trunc-p
      (format stream "~&..."))))

(defmacro define-matrix-like-printer ((instance class stream) slots &body body)
  "Simple wrapper for the standard printing style."
  (check-type class symbol)
  (check-type instance symbol)
  (check-type stream symbol)
  `(defmethod print-object ((,instance ,class) ,stream)
     (print-unreadable-object (,instance ,stream :type t)
       (format ,stream "~s with ~a x ~a elements~&" (lla-type ,instance) (nrow ,instance) (ncol ,instance))
       (with-slots ,slots ,instance
         ,@body))))

(define-matrix-like-printer (matrix dense-matrix-like stream) ()
  (print-matrix matrix stream))

(define-matrix-like-printer (matrix upper-triangular-matrix stream) ()
  (print-matrix matrix stream :mask #'upper-triangular-mask))

(define-matrix-like-printer (matrix lower-triangular-matrix stream) ()
  (print-matrix matrix stream :mask #'lower-triangular-mask))

(define-matrix-like-printer (matrix hermitian-matrix stream) ()
  (print-matrix matrix stream :mask #'hermitian-mask))
