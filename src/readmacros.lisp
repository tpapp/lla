(in-package :lla)

(defun read-string-upto-char (char &optional (stream *standard-input*)
                                   (eof-error-p t) eof-value recursive-p)
  "Read from stream until encountering CHAR.  Return everything before
as a string."
  (handler-case
      (with-output-to-string (out)
        (loop for c = (read-char stream t nil recursive-p)
              until (char= c char)
              do (write-char c out)))
    (end-of-file (e)
      (if eof-error-p
          (error e)
          eof-value))))

(defun read-vector-or-matrix (stream c n)
  "Function for read macro to read LLA vectors or matrices.  Syntax:
#[n]v[type][:kind](e1 e2 ...), where type can be i (integer), s or
d (single/double), cs or cd (complex single/double).  If n is present,
elements are interpeted as a matrix with n columns (0 gives a row
matrix).  Kind gives the matrix kind for matrices (dense, upper,
lower, hermitian).  Examples:

 #v(1 2 3)               vector of 3 elements
 #2v(1 2 3 4)            2x2 matrix
 #0vcs(1 2 3)            1x3 matrix of complex single elements
 #2v:hermitian(1 2 3 4)  2x2 hermitian matrix (only upper triangle is used).
 #v:diagonal(1 2 3 4)    4x4 diagonal matrix."
  (declare (ignore c)
           (optimize debug))
  (bind ((specification (read-string-upto-char #\( stream t nil t))
         (elements (read-delimited-list #\) stream t))
         ((:values type-string kind-string) 
          (aif (position #\: specification)
               (values (subseq specification 0 it) (subseq specification (1+ it)))
               (values specification nil)))
         (type (alexandria:eswitch (type-string :test #'string-equal)
                 ("" nil)
                 ("s" :single)
                 ("d" :double)
                 ("cs" :complex-single)
                 ("cd" :complex-double)
                 ("i" :integer))))
    (cond
      ((null elements) (error "no elements provided"))
      (n (create-matrix n elements :lla-type type :kind
                        (if kind-string
                            (alexandria:eswitch (kind-string :test #'string-equal)
                              ("dense" :dense)
                              ("upper" :upper-triangular)
                              ("lower" :lower-triangular)
                              ("hermitian" :hermitian))
                            :dense)))
      (kind-string (if (string-equal kind-string "diagonal")
                       (create-diagonal elements type)
                       (error "unrecognized matrix/vector kind")))
      (t (create-nv elements type)))))

(defreadtable lla:v-syntax
  (:merge :standard)
  (:dispatch-macro-char #\# #\v #'read-vector-or-matrix))
