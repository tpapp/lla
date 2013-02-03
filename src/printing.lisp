(in-package :lla)

;;;; Printing and formatting

;;;; General variables and utility functions

(defun print-length-truncate (dimension)
  "Return values (min dimension *print-length*) and whether the
constraint is binding."
  (if (or (not *print-length*) (<= dimension *print-length*))
      (values dimension nil)
      (values *print-length* t)))

(defvar *lla-print-precision* 5
  "number of digits after the decimal point when printing numeric matrices")

(defun standard-numeric-formatter (x)
  "Standard formatter for matrix printing.  Respects
*print-lla-precision*, and formats complex numbers as a+bi, eg
0.0+1.0i."
  ;; ?? do we want a complex numbers to be aligned on the +, like R? I
  ;; am not sure I like that very much, and for a lot of data, I would
  ;; visualize it graphically anyhow (I hate tables of 7+ numbers in
  ;; general).  -- Tamas, 2009-sep-13
  (let ((precision *lla-print-precision*))
    (typecase x
      (integer (format nil "~d" x))
      (real (format nil "~,vf" precision x))
      (complex (format nil "~,vf+~,vfi"
                       precision (realpart x)
                       precision (imagpart x)))
      (t (format nil "~a" x)))))

;;;; matrices

;;; None of the code below assumes anything about classes, it just
;;; uses the mref interface to access elements.

(defvar *lla-print-matrix-aligned* t
  "If non-nil, characters will be aligned.")

(defvar *lla-print-matrix-paddig* 1
  "Number of spaces between columns.")

(defun print-matrix (matrix stream masked-element
                     &key (formatter #'standard-numeric-formatter))
  "Format and print the elements of matrix to stream, using
*LLA-PRINT-MATRIX-PADDING* spaces between columns.  If
*LLA-PRINT-MATRIX-ALIGNED*, columns will be right-aligned.  Prints at most
*PRINT-LENGTH* rows and columns, indicating more with a ...  Uses MREF for
element access, printing MASKED-ELEMENT for masked elements.."
  ;; ?? maybe column & row labels, not a high priority at the moment
  (let+ (((&values nrow row-trunc?) (print-length-truncate (aops:nrow matrix)))
	 ((&values ncol col-trunc?) (print-length-truncate (aops:ncol matrix)))
	 (formatted-elements (make-array (list nrow ncol)))
	 (column-widths (make-array ncol :element-type 'fixnum
                                         :initial-element 0))
	 (padding (make-array *lla-print-matrix-paddig*
                              :element-type 'character
                              :initial-element #\space))
	 (aligned? *lla-print-matrix-aligned*))
    ;; first pass - format elements, measure width
    (dotimes (col ncol)
      (dotimes (row nrow)
	(let+ (((&values element masked?) (mref matrix row col))
               (formatted-element (if masked?
                                      masked-element
                                      (funcall formatter element)))
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
	  (if aligned?
	      (format stream "~V@A" (aref column-widths col) elt)
	      (princ elt stream))))
      (when col-trunc?
	(princ padding stream)
	(princ "..." stream)))
    (when row-trunc?
      (format stream "~&..."))))
