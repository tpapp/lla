;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla)

;;; Helper functions


(eval-when (:compile-toplevel :load-toplevel :execute)

  (defun value-to-memory% (lla-type)
    "Return a (LAMBDA (POINTER INDEX VALUE) ...) form that can be used to
write an element to an array in memory."
    `(lambda (pointer index value)
       ,(ecase lla-type
          (:single
           `(setf (mem-aref pointer :float index)
                  (coerce value 'single-float)))
          (:double
           `(setf (mem-aref pointer :double index)
                  (coerce value 'double-float)))
          (:complex-single
           `(let ((index2 (the fixnum (* 2 index))))
              (setf (mem-aref pointer :float index2)
                    (coerce (realpart value) 'single-float)
                    (mem-aref pointer :float (the fixnum (1+ index2)))
                    (coerce (imagpart value) 'single-float))))
          (:complex-double
           `(let ((index2 (the fixnum (* 2 index))))
              (setf (mem-aref pointer :double index2)
                    (coerce (realpart value) 'double-float)
                    (mem-aref pointer :double (the fixnum (1+ index2)))
                    (coerce (imagpart value) 'double-float))))
          (:integer
           `(setf (mem-aref pointer #-lla::int64 :int32
                                    #+lla::int64 :int64
                                    index)
                  (coerce value
                          '(signed-byte #-lla::int64 32 #+lla::int64 64)))))
       (values)))
  
  (defun value-from-memory% (lla-type)
    "Return a (LAMBDA (POINTER INDEX) ...) form that can be used to read an
element from an array in memory."
    `(lambda (pointer index)
       ,(ecase lla-type
          (:single
           `(mem-aref pointer :float index))
          (:double
           `(mem-aref pointer :double index))
          (:complex-single
           `(let ((index2 (the fixnum (* 2 index))))
              (complex (mem-aref pointer :float index2)
                       (mem-aref pointer :float (1+ index2)))))
          (:complex-double
           `(let ((index2 (the fixnum (* 2 index))))
              (complex (mem-aref pointer :double index2)
                       (mem-aref pointer :double (1+ index2)))))
          (:integer
           `(mem-aref pointer #-lla::int64 :int32
                              #+lla::int64 :int64
                              index)))))
  
  (defun expand-specifications% (clause specifications)
    "Expand specifications using (clause lla-type element-type)."
    (mapcan (lambda+ ((lla-type &rest element-types))
              (mapcar (curry clause lla-type) element-types))
            specifications))

  (defun all-to-specifications% ()
    '((:single single-float *)
      (:double single-float double-float *)
      (:complex-single single-float double-float
       (complex single-float) *)
      (:complex-double single-float double-float
       (complex single-float)
       (complex double-float) *)
      (:integer (signed-byte #-lla::int64 32 #+lla::int64 64) *)))

  (defun all-from-specifications% ()
    '((:single lla-single *)
      (:double lla-double *)
      (:complex-single lla-complex-single *)
      (:complex-double lla-complex-double *)
      (:integer lla-integer *))))

(eval-when (:compile-toplevel :load-toplevel :execute)
  (defmacro copy-to-memory% (array pointer lla-type
                             &optional (specifications (all-to-specifications%)))
    "Helper function for generating optimized code conditioning on lla-type and
the type of array."
    (let+ (((&once-only array pointer lla-type))
           ((&with-gensyms index size))
           ((&flet clause (lla-type% element-type%)
              (let* ((simple? (not (eq element-type% '*)))
                     (array-type% `(,(if simple? 'simple-array 'array)
                                     ,element-type% *)))
                `((and (typep ,array ',array-type%)
                       (eq ,lla-type ,lla-type%))
                  (locally
                      (declare ,@(when simple? '((optimize speed (safety 0))))
                               (type ,array-type% ,array))
                    (loop for ,index below ,size do
                      (,(value-to-memory% lla-type%) ,pointer ,index
                       (row-major-aref ,array ,index)))))))))
      `(let ((,size (array-total-size ,array)))
         (cond
           ,@(expand-specifications% #'clause specifications)
           (t (error "Don't know how to perform operation.")))))))

(defmacro copy-from-memory% (array pointer lla-type
                             &optional (specifications (all-from-specifications%)))
  "Helper function for generating optimized code conditioning on lla-type and
the type of array."
  (let+ (((&once-only array pointer lla-type))
         ((&with-gensyms index size))
         ((&flet clause (lla-type% element-type%)
            (let* ((simple? (not (eq element-type% '*)))
                   (array-type% `(,(if simple? 'simple-array 'array)
                                   ,element-type% *)))
              `((and (typep ,array ',array-type%)
                     (eq ,lla-type ,lla-type%))
                (locally
                    (declare ,@(when simple? '((optimize speed (safety 0))))
                             (type ,array-type% ,array))
                  (loop for ,index below ,size do
                    (setf (row-major-aref ,array ,index)
                          (coerce
                           (,(value-from-memory% lla-type%)
                             ,pointer ,index)
                           ',(if (eq element-type% '*)
                                 t
                                 element-type%))))))))))
    `(let ((,size (array-total-size ,array)))
       (cond
         ,@(expand-specifications% #'clause specifications)
         (t (error "Don't know how to perform operation."))))))

(defun copy-array-to-memory (array pointer lla-type)
  "Copy the contents of ARRAY to the memory area of type LLA-TYPE at POINTER."
  (check-type array array)
  (copy-to-memory% array pointer lla-type)
  (values))

(defun copy-array-from-memory (array pointer lla-type)
  "Copy the memory area of type LLA-TYPE at POINTER to ARRAY."
  (check-type array array)
  (copy-from-memory% array pointer lla-type)
  (values))

(defun create-array-from-memory (pointer lla-type dimensions
                                 &optional (element-type (lla-to-lisp-type
                                                          lla-type)))
  "Create an array from contents at POINTER."
  (aprog1 (make-array dimensions :element-type element-type)
    (copy-array-from-memory it pointer lla-type)))

(defmacro transpose-to-memory% (matrix pointer lla-type
                                &optional (specifications (all-to-specifications%)))
  "Helper function for generating optimized code conditioning on lla-type and
the type of the matrix."
  (let+ (((&once-only matrix pointer lla-type))
         ((&with-gensyms index nrow ncol row-index col-index))
         ((&flet clause (lla-type% element-type%)
            (let* ((simple? (not (eq element-type% t)))
                   (array-type% `(,(if simple? 'simple-array 'array)
                                   ,element-type% (* *))))
              `((and (typep ,matrix ',array-type%)
                     (eq ,lla-type ,lla-type%))
                (locally (declare ,@(when simple?
                                      '((optimize speed (safety 0))))
                                  (type ,array-type% ,matrix))
                  (loop for ,col-index fixnum below ,ncol do
                    (loop for ,row-index fixnum below ,nrow do
                      (,(value-to-memory% lla-type%) ,pointer ,index
                       (aref ,matrix ,row-index ,col-index))
                      (incf ,index)))))))))
    `(let+ (((,nrow ,ncol) (array-dimensions ,matrix))
            (,index 0))
       (declare (type fixnum ,nrow ,ncol ,index))
       (cond
         ,@(expand-specifications% #'clause specifications)
         (t (error "Don't know how to perform operation."))))))

(defmacro transpose-from-memory% (matrix pointer lla-type
                                  &optional (specifications
                                             (all-from-specifications%)))
  "Helper function for generating optimized code conditioning on lla-type and
the type of the matrix."
  (let+ (((&once-only matrix pointer lla-type))
         ((&with-gensyms index nrow ncol row-index col-index))
         ((&flet clause (lla-type% element-type%)
            (let* ((simple? (not (eq element-type% t)))
                   (array-type% `(,(if simple? 'simple-array 'array)
                                   ,element-type% (* *))))
              `((and (typep ,matrix ',array-type%)
                     (eq ,lla-type ,lla-type%))
                (locally (declare ,@(when simple?
                                      '((optimize speed (safety 0))))
                                  (type ,array-type% ,matrix))
                  (loop for ,col-index fixnum below ,ncol do
                    (loop for ,row-index fixnum below ,nrow do
                      (setf (aref ,matrix ,row-index ,col-index)
                            (coerce
                             (,(value-from-memory% lla-type%) ,pointer ,index)
                             ',(if (eq element-type% '*)
                                   t
                                   element-type%)))
                      (incf ,index)))))))))
    `(let+ (((,nrow ,ncol) (array-dimensions ,matrix))
            (,index 0))
       (declare (type fixnum ,nrow ,ncol ,index))
       (cond
         ,@(expand-specifications% #'clause specifications)
         (t (error "Don't know how to perform operation."))))))

(defun transpose-matrix-to-memory (matrix pointer lla-type)
  "Transpose the contents of ARRAY to the memory area of type LLA-TYPE at
  POINTER."
  (check-type matrix (array * (* *)))
  (transpose-to-memory% matrix pointer lla-type))

(defun transpose-matrix-from-memory (matrix pointer lla-type)
  "Transpose the contents of ARRAY to the memory area of type LLA-TYPE at
  POINTER."
  (check-type matrix (array * (* *)))
  (transpose-from-memory% matrix pointer lla-type))

(defun create-transposed-matrix-from-memory (pointer lla-type dimensions
                                             &optional (element-type
                                                        (lla-to-lisp-type
                                                         lla-type)))
  "Create a matrix from transposed contents at POINTER."
  (aprog1 (make-array dimensions :element-type element-type)
    (transpose-matrix-from-memory it pointer lla-type)))

