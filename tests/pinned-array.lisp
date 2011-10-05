;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla-tests)

(deftestsuite pinned-array-tests (lla-tests)
  ())

;;;; Basic testing of copying functions
;;; 
;;; Should not require LLA facilities other than those tested.

(defun flatten-complex-vector (vector)
  "Separate real and complex elements in a vector and stack them
consecutively."
  (coerce (iter
            (for element :in-vector vector)
            (collecting (realpart element))
            (collecting (imagpart element)))
          'vector))

(defun check-memory-contents (pointer cffi-type values)
  "Check that the contents of the memory (of CFFI-TYPE) at POINTER is the same
as VALUES.  VALUES should be flattened if complex."
  (loop
    for index from 0
    for value across values
    always (let* ((value-in-memory (mem-aref pointer cffi-type index))
                  (same? (= value value-in-memory)))
             (unless same?
               (warn "value (type ~A) at index ~A is ~A, expected ~A."
                     cffi-type index value-in-memory value))
             same?)))

(defun check-memory-contents2 (pointer internal-type values)
  "Like CHECK-MEMORY-CONTENTS, with autoconversion of values and calculation
of INTERNAL-TYPE."
  (check-memory-contents pointer
                         (alexandria:eswitch (internal-type)
                           (lla::+single+ :float)
                           (lla::+double+ :double)
                           (lla::+complex-single+ :float)
                           (lla::+complex-double+ :double)
                           (lla::+integer+ #-lla::int64 :int32
                                           #+lla::int64 :int64))
                         (let ((vector (aetypecase values
                                         (vector it)
                                         (array (flatten-array it))
                                         (t (coerce it 'vector)))))
                           (if (lla::complex? internal-type)
                               (flatten-complex-vector vector)
                               vector))))

(addtest (pinned-array-tests)
  copy-to-memory
  (with-foreign-temporary-buffer (pointer 100)
    (let+ (((&flet check-copy (array internal-type)
              (lla::copy-array-to-memory array pointer internal-type)
              (check-memory-contents2 pointer internal-type array)))
           ((&macrolet ensure-copy (element-type internal-type)
              `(ensure (check-copy (random-array ,element-type 5 7)
                                   ,internal-type))))
           ((&macrolet ensure-copy-error (element-type internal-type)
              `(ensure-error
                 (lla::copy-array-to-memory (random-array ,element-type 19)
                                            pointer ,internal-type)))))
      ;; single float
      (ensure-copy t lla::+single+)
      (ensure-copy 'lla-integer lla::+single+)
      (ensure-copy 'lla-single lla::+single+)
      (ensure-copy 'lla-double lla::+single+)
      (ensure-copy-error 'lla-complex-single lla::+single+)
      ;; double float
      (ensure-copy t lla::+double+)
      (ensure-copy 'lla-integer lla::+double+)
      (ensure-copy 'lla-single lla::+double+)
      (ensure-copy 'lla-double lla::+double+)
      (ensure-copy-error 'lla-complex-single lla::+double+)
      ;; complex single
      (ensure-copy t lla::+complex-single+)
      (ensure-copy 'lla-integer lla::+complex-single+)
      (ensure-copy 'lla-single lla::+complex-single+)
      (ensure-copy 'lla-double lla::+complex-single+)
      (ensure-copy 'lla-complex-single lla::+complex-single+)
      (ensure-copy 'lla-complex-double lla::+complex-single+)
      ;; complex double
      (ensure-copy t lla::+complex-double+)
      (ensure-copy 'lla-integer lla::+complex-double+)
      (ensure-copy 'lla-single lla::+complex-double+)
      (ensure-copy 'lla-double lla::+complex-double+)
      (ensure-copy 'lla-complex-single lla::+complex-double+)
      (ensure-copy 'lla-complex-double lla::+complex-double+))))

(addtest (pinned-array-tests)
  transpose-to-memory
  (with-foreign-temporary-buffer (pointer 100)
    (let+ (((&flet check-transpose (matrix internal-type)
              (lla::transpose-matrix-to-memory matrix pointer internal-type)
              (check-memory-contents2 pointer internal-type (transpose matrix))))
           ((&macrolet ensure-transpose (element-type internal-type)
              `(ensure (check-transpose (random-array ,element-type 5 6)
                                        ,internal-type))))
           ((&macrolet ensure-transpose-error (element-type internal-type)
              `(ensure-error 
                 (lla::transpose-matrix-to-memory 
                  (random-array ,element-type 7 5) pointer ,internal-type)))))
      ;; single float
      (ensure-transpose t lla::+single+)
      (ensure-transpose 'lla-integer lla::+single+)
      (ensure-transpose 'lla-single lla::+single+)
      (ensure-transpose 'lla-double lla::+single+)
      (ensure-transpose-error 'lla-complex-single lla::+single+)
      ;; ;; double float
      (ensure-transpose t lla::+double+)
      (ensure-transpose 'lla-integer lla::+double+)
      (ensure-transpose 'lla-single lla::+double+)
      (ensure-transpose 'lla-double lla::+double+)
      (ensure-transpose-error 'lla-complex-single lla::+double+)
      ;; complex single
      (ensure-transpose t lla::+complex-single+)
      (ensure-transpose 'lla-integer lla::+complex-single+)
      (ensure-transpose 'lla-single lla::+complex-single+)
      (ensure-transpose 'lla-double lla::+complex-single+)
      (ensure-transpose 'lla-complex-single lla::+complex-single+)
      (ensure-transpose 'lla-complex-double lla::+complex-single+)
      ;; complex double
      (ensure-transpose t lla::+complex-double+)
      (ensure-transpose 'lla-integer lla::+complex-double+)
      (ensure-transpose 'lla-single lla::+complex-double+)
      (ensure-transpose 'lla-double lla::+complex-double+)
      (ensure-transpose 'lla-complex-single lla::+complex-double+)
      (ensure-transpose 'lla-complex-double lla::+complex-double+))))

(addtest (pinned-array-tests)
  copy-from-memory
  (with-foreign-temporary-buffer (pointer 100)
    (let+ (((&flet check-copy-to-and-from (element-type internal-type)
              (let ((array (random-array element-type 5 8)))
                (lla::copy-array-to-memory array pointer internal-type)
                (equalp array (lla::create-array-from-memory
                               pointer internal-type (array-dimensions array))))))
           ((&macrolet ensure-to-and-from (element-type internal-type)
              `(ensure (check-copy-to-and-from ,element-type ,internal-type)))))
      (ensure-to-and-from 'lla-single lla::+single+)
      (ensure-to-and-from 'lla-double lla::+double+)
      (ensure-to-and-from 'lla-complex-single lla::+complex-single+)
      (ensure-to-and-from 'lla-complex-double lla::+complex-double+)
      (ensure-to-and-from 'lla-integer lla::+integer+))))

(addtest (pinned-array-tests)
  transpose-from-memory
  (with-foreign-temporary-buffer (pointer 100)
    (let+ (((&flet check-transpose-to-and-from (element-type internal-type)
              (let ((matrix (random-array element-type 5 8)))
                (lla::transpose-matrix-to-memory matrix pointer internal-type)
                (equalp matrix (lla::create-transposed-matrix-from-memory
                                pointer internal-type (array-dimensions matrix))))))
           ((&macrolet ensure-to-and-from (element-type internal-type)
              `(ensure (check-transpose-to-and-from ,element-type ,internal-type)))))
      (ensure-to-and-from 'lla-single lla::+single+)
      (ensure-to-and-from 'lla-double lla::+double+)
      (ensure-to-and-from 'lla-complex-single lla::+complex-single+)
      (ensure-to-and-from 'lla-complex-double lla::+complex-double+)
      (ensure-to-and-from 'lla-integer lla::+integer+))))

;;;; Pinning tests
;;; 
;;; Test the pinning framework, regardless of its implementation (copying,
;;; sharing, etc).

(defun coercible? (source-type target-type)
  "Permitted coercions for internal LLA types.  It is guaranteed that
CL:COERCE can perform these coercions on the corresponding lisp types. Basic
summary: (1) integers can be coerced to anything, (2) single<->double
precision coercions are possible both ways, (3) real floats can be upgraded to
complex."
  (cond
    ;; always valid
    ((eq source-type lla::+integer+) t)
    ;; nothing else can be converted to integer
    ((eq target-type lla::+integer+) nil)
    ;; no complex->real
    ((and (lla::complex? source-type)
          (not (lla::complex? target-type)))
     nil)
    ;; all the rest should be possible
    (t t)))

(defun coercible-pairs-list ()
  ;;   "Generate the list of all LLA (source target) pairs for which
  ;; coercible? holds.  For internal use only, NOT EXPORTED."
  (let (coercible)
    (dolist (source lla::+internal-types+)
      (dolist (target lla::+internal-types+)
        (when (coercible? source target)
          (push (list source target) coercible))))
    ;; reverse only for cosmetic purposes (debugging)
    (nreverse coercible)))

(defun test-pinning-readonly (source-type destination-type
                              &key (length 50))
  "Test array pinning (input only)."
  (let ((vector (random-array (lla::lisp-type source-type) length)))
    (lla::with-pinned-array (pointer vector destination-type nil nil nil nil)
      (check-memory-contents2 pointer destination-type vector))))

(defun test-pinning-copy (source-type destination-type
                          &key (length 50))
  "Test array pinning (copy requested).  The memory is modified, compared
again, and the original vector is checked to ensure that it remains constant."
  (let* ((vector (random-array (lla::lisp-type source-type) length))
         (copy (copy-seq vector))
         (copy-inc (map 'vector #'1+ copy)))
    (lla::with-pinned-array (pointer vector destination-type nil
                                     :copy nil nil)
      ;; check equality
      (unless (check-memory-contents2 pointer destination-type copy)
        (return-from test-pinning-copy nil))
      ;; increment memory area, check equality, and that original is intact
      (dotimes (i length)
        (incf (lla::foreign-aref pointer destination-type i) 1))
      (and (check-memory-contents2 pointer destination-type copy-inc)
           (equalp vector copy)))))

(defun test-pinning-output (source-type destination-type
                            &key (length 50))
  "Test array pinning (with output).  Procedure is similar to
test-pinning-copy, except that the output is also checked."
  (let* ((vector (random-array (lla::lisp-type source-type) length))
         vector-inc
         (copy (copy-seq vector))
         (copy-inc (map 'vector #'1+ copy)))
    (lla::with-pinned-array (pointer vector destination-type nil
                                     vector-inc nil nil)
      ;; check equality
      (unless (check-memory-contents2 pointer destination-type copy)
        (return-from test-pinning-output nil))
      ;; increment memory area, check equality, and that original is intact
      (dotimes (i length)
        (incf (lla::foreign-aref pointer destination-type i) 1))
      (unless (and (check-memory-contents2 pointer destination-type copy-inc)
                   (equalp vector copy))
        (return-from test-pinning-output nil)))
    ;; check output
    (every #'= vector-inc copy-inc)))

(defun test-vector-output (internal-type &key (length 50))
  (let ((vector (random-array (lla::lisp-type internal-type) length))
        output)
    (lla::with-array-output (pointer output internal-type length nil)
      (dotimes (index length)
        (setf (lla::foreign-aref pointer internal-type index)
              (aref vector index))))
    (equalp vector output)))


;; pinning, input only

(addtest (pinned-array-tests)
  pinning-input
  (iter 
    (for (source destination) :in (coercible-pairs-list))
    (ensure
     (test-pinning-readonly source destination))))

;; pinning, copy

(addtest (pinned-array-tests)
  pinning-copy
  (iter 
    (for (source destination) :in (coercible-pairs-list))
    (ensure
     (test-pinning-copy source destination))))

;; pinning, output

(addtest (pinned-array-tests)
  pinning-output
  (iter 
    (for (source destination) :in (coercible-pairs-list))
    (ensure
     (test-pinning-output source destination))))

;; vector-output

(addtest (pinned-array-tests)
  vector-output
  (ensure (test-vector-output lla::+integer+))
  (ensure (test-vector-output lla::+single+))
  (ensure (test-vector-output lla::+double+))
  (ensure (test-vector-output lla::+complex-single+))
  (ensure (test-vector-output lla::+complex-double+)))
