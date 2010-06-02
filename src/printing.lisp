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

;;;; diagonals

(defmethod print-object ((obj diagonal) stream)
  (print-unreadable-object (obj stream :type t)
    (with-slots (elements) obj
      (format stream "~s with ~a elements: ~a" (array-lla-type elements)
              (length elements)
              elements))))

;;;; matrices

;;; None of the code below assumes anything about classes, it just
;;; uses the mref interface to access elements.

(defvar *print-matrix-aligned* t
  "If non-nil, characters will be aligned.")

(defvar *print-matrix-paddig* 1
  "Number of spaces between columns.")

(defgeneric matrix-masked-element-string (kind)
  (:documentation "Return a string that can be used to represend
  masked elements when printing.  Return value may not be modified.")
  (:method ((kind (eql :dense)))
    nil)
  (:method ((kind (eql :upper)))
    ".")
  (:method ((kind (eql :lower)))
    ".")
  (:method ((kind (eql :hermitian)))
    "*"))

(defun print-matrix (matrix stream &key 
		     (formatter #'standard-numeric-formatter)
                     (accessor #'mref)
                     (mask (constantly nil)))
  "Format and print the elements of matrix to stream, using
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
      (format stream "  ")
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

(defmethod print-object ((matrix dense-matrix-like) stream)
  (bind (((:slots-r/o elements nrow ncol) matrix)
         (kind (matrix-kind matrix))
         (masked-element-string (matrix-masked-element-string kind)))
    (print-unreadable-object (matrix stream :type t)
      (format stream "~s with ~a x ~a elements~&"
              (array-lla-type elements) nrow ncol)
      (print-matrix matrix stream 
                    :mask (lambda (row col)
                            (if (matrix-mask kind row col)
                                nil
                                masked-element-string))))))
