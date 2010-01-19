(in-package :lla)

(defun read-vector-or-matrix (stream c n)
  "Function for read macro to read LLA vectors or matrices.  Syntax:
#[n]v[type][:kind](e1 e2 ...), where type can be i (integer), s or
d (single/double), cs or cd (complex single/double).  If n is present,
elements are interpeted as a matrix with n columns (0 gives a row
matrix).  Kind gives the matrix kind for matrices (dense,
upper-triangular, lower-triangular, hermitian).  Examples

 #v(1 2 3)               vector of 3 elements
 #2v(1 2 3 4)            2x2 matrix
 #0vcs(1 2 3)            1x3 matrix of complex single elements
 #2v:hermitian(1 2 3 4)  2x2 hermitian matrix (only upper triangle is used)."
  (declare (ignore c))
  (let (type-string
        kind-string
        collected-chars
        (state :init))
    (macrolet ((save-collected-to (place)
                 "Save collected characters (if any) to place."
                 `(when collected-chars
                    (setf ,place (coerce (nreverse collected-chars) 'string)))))
      (loop
        (let ((char (read-char stream t nil t)))
          (case char
            (#\: (ecase state           ; #vtype:
                   (:init (save-collected-to type-string)
                      (setf collected-chars nil
                            state :expecting-kind))
                   (:expecting-kind (error "Two :'s in specification."))))
            (#\( (ecase state
                   (:init (save-collected-to type-string)) ; #vtype(
                   (:expecting-kind (save-collected-to kind-string))) ; #vtype:kind(
               (return))
            (otherwise (push char collected-chars))))))
    (flet ((find-or-default (string-or-nil default argument-description alist)
             "Identify string when non-nil using ALIST, otherwise
             return DEFAULT.  ARGUMENT-DESCRIPTION is for the error message."
             (if string-or-nil
                 (let ((entry (assoc string-or-nil alist :test #'string=)))
                   (if entry
                       (cdr entry)
                       (error "Could not identify ~A as ~S argument."
                              string-or-nil argument-description)))
                 default)))
      (let ((lla-type (find-or-default type-string nil "an element type"
                                       '(("i" . :integer)
                                         ("s" . :single)
                                         ("d" . :double)
                                         ("cs" . :complex-single)
                                         ("cd" . :complex-double))))
            (elements (read-delimited-list #\) stream t)))
        (cond
          (n (create-matrix n elements :lla-type lla-type :kind
                            (find-or-default kind-string :dense "a matrix kind"
                                             '(("dense" . :dense)
                                               ("upper-triangular" . :upper-triangular)
                                               ("lower-triangular" . :lower-triangular)
                                               ("hermitian" . :hermitian)))))
          ((string= kind-string "diagonal") (create-diagonal elements lla-type))
          (kind-string (error "vectors don't have kinds."))
          (t (create-nv elements lla-type)))))))

(defreadtable lla-readtable
  (:merge :standard)
  (:dispatch-macro-char #\# #\v #'read-vector-or-matrix))
