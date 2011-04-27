(in-package :lla)

(define-make-symbol% :lla)

(deftype symbol* () '(and symbol (not null)))
(defun symbolp* (object) (typep object 'symbol*))

(deftype matrix (&optional (m '*) (n '*))
  `(array * (,m ,n)))

(deftype matrix? (object)
  "Return a boolean indicating whether OBJECT is a MATRIX."
  (typep object 'matrix))

(defun square? (matrix)
  "Test if a matrix (in the generalized sense, ie an object that has nrow and ncol)
is square."
  (= (nrow matrix) (ncol matrix)))

(defmacro row-major-loop ((matrix index row col &key
                                 (nrow (gensym* '#:nrow)) (ncol (gensym* '#:ncol)))
                          &body body)
  "Loop through row-major MATRIX, incrementing INDEX (flat row-major index), ROW, and
COL."
  (check-types (row col index nrow ncol) symbol)
  (once-only (matrix)
    `(bind (((,nrow ,ncol) (array-dimensions ,matrix))
            (,index 0))
       (dotimes (,row ,nrow)
         (dotimes (,col ,ncol)
           ,@body
           (incf ,index))))))

(defmacro define-ondemand-slot ((instance-and-class slot-name) &body body)
  `(defmethod slot-missing (,(gensym) ,instance-and-class (slot-name (eql ',slot-name))
                            (operation (eql 'slot-value)) &optional new-value)
     (declare (ignore new-value))
     ,@body))

(defun transposed-dimensions (m n transpose &optional (library :blas))
  "Process dimensions of a matrix, which may be transposed.  Third value is the enum
expected by the library."
  (ecase transpose
    ((nil) (values m n (ecase library (:blas :CblasNoTrans) (:lapack (char-code #\N)))))
    ((t) (values n m (ecase library (:blas :CblasTrans) (:lapack #\T))))
    ((*) (values n m (ecase library (:blas :CblasConjTrans) (:lapack #\C))))))

(defun maybe-vector-as-matrix (vector-or-matrix orientation &optional 
                               transpose (library :blas))
  "Process dimensions of vectors which are meant to be reoriented as matrices.
Orientation is :ROW or :COLUMN, and TRANSPOSE indicates transposition (see TRANS).
Return (values DIMENSION0 DIMENSION1 ORIENTATION LEADING-DIMENSION TRANSPOSE-ENUM).
ORIENTATION is NIL or the ORIENTATION argument to the function, depending on whether
the argument was a vector."
  (bind (((d1 &optional d2 d-rest) (array-dimensions vector-or-matrix))
         ((:flet error% ()) (error "~A is not a vector or a matrix." vector-or-matrix))
         ((:values m n orientation)
          (cond
            (d-rest (error%))
            (d2 (values d1 d2))
            (d1 (ecase orientation
                  (:row (values 1 d1 :row))
                  (:column (values d1 1 :column))))
            (t (error%))))
         ((:values m n transpose-enum) (transposed-dimensions m n transpose library)))
    (values m n orientation n transpose-enum)))

(defun vector-or-matrix-dimensions (m n orientation)
  "Return dimensions for matrices that can potentially be vectors, depending on
 orientation."
  (ecase orientation
    (:row (assert (= m 1)) n)
    (:column (assert (= n 1)) m)
    ((nil) (list m n))))

(defun maybe-pick-first-element (array pick?)
  "When PICK"
  (if pick?
      (progn (assert (= (array-total-size array) 1) ()
                     "~A is supposed to have only one element." array)
             (row-major-aref array 0))
      array))

(defun matrix-from-first-rows (matrix nrow orientation)
  "Create a matrix (or vector, depending on ORIENTATION) from the first rows NRHS of
MATRIX.  Used for interfacing with xGELS, extracting R from QR decompositions, etc."
  (bind (((nil n) (array-dimensions matrix)))
    (copy-array (displace-array matrix 
                                (vector-or-matrix-dimensions nrow n orientation)))))

;;;; this is missing from CFFI at the moment
(define-with-multiple-bindings with-foreign-pointer)

;; ;; #+sbcl (eval-when (:compile-toplevel :load-toplevel :execute)
;; ;;          (pushnew :muffle-notes cl:*features*))

;; (declaim (inline as-scalar%))
;; (defun as-scalar% (vector)
;;   "Pick the first element from a vector.  No checking."
;;   (row-major-aref vector 0))

;; ;;; muffling notes
;; (defmacro muffle-optimization-notes (&body body)
;;   "This macro silences compiler optimization notes."
;;   `(locally 
;;        #+sbcl (declare (sb-ext:muffle-conditions sb-ext:compiler-note))
;;        ,@body))

;; (defun maybe-wrap-list (prefix body)
;;   "When PREFIX, wraps BODY in a list, starting with PREFIX, otherwise
;; return BODY.  Intended for use in macros."
;;   (if prefix
;;       (list (concatenate 'list prefix body))
;;       body))

;; (defun maybe-wrap (prefix &rest body)
;;   "When PREFIX, wraps BODY in a list, starting with PREFIX, otherwise
;; return BODY.  Intended for use in macros."
;;   (maybe-wrap-list prefix body))

;; ;;; array utilities

;; (declaim (inline zero-like))
;; (defun zero-like (array)
;;   "Return 0 coerced to the element type of ARRAY."
;;   (coerce 0 (array-element-type array)))

;; ;;; unfortunately, simple-vector is already taken by CL, for element
;; ;;; types T, so we define the new type SIMPLE-ARRAY1
;; (deftype simple-array1 (&optional element-type length)
;;   `(simple-array ,element-type (,length)))

;; (defun simple-array? (object)
;;   (typep object 'simple-array))

;; (defun simple-array1? (object)
;;   (typep object 'simple-array1))

;; (defun as-simple-array1 (array)
;;   "Return elements of ARRAY as a SIMPLE-ARRAY1.  Array is not
;; necessarily copied if it is already of the correct type."
;;   (etypecase array
;;     (simple-array1 array)
;;     (array (copy-seq (displace-array array (array-total-size array))))))
