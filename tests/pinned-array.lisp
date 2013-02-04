;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla-tests)

(defsuite pinned-array-suite (tests))

;;;; Basic testing of copying functions
;;;
;;; Should not require LLA facilities other than those tested.

(defun flatten-complex-vector (vector)
  "Separate real and complex elements in a vector and stack them consecutively."
  (coerce (loop
            for element across vector
            collect (realpart element)
            collect (imagpart element))
          'vector))

(defun check-memory-contents (pointer cffi-type values)
  "Check that the contents of the memory (of CFFI-TYPE) at POINTER is the same as VALUES.  VALUES should be flattened if complex."
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
  "Like CHECK-MEMORY-CONTENTS, with autoconversion of values and calculation of INTERNAL-TYPE."
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
                                         (array (aops:flatten it))
                                         (t (coerce it 'vector)))))
                           (if (lla::complex? internal-type)
                               (flatten-complex-vector vector)
                               vector))))

(deftest copy-to-memory (pinned-array-suite)
  (with-foreign-temporary-buffer (pointer 100)
    (let+ (((&flet check-copy (array internal-type)
              (lla::copy-array-to-memory array pointer internal-type)
              (assert-true (check-memory-contents2 pointer internal-type array))))
           ((&flet assert-same-copy (element-type internal-type)
              (check-copy (random-array element-type 5 7)
                          internal-type)))
           ((&flet assert-copy-error (element-type internal-type)
              (assert-condition error
                  (lla::copy-array-to-memory (random-array element-type 19)
                                             pointer internal-type)))))
      ;; single float
      (assert-same-copy t lla::+single+)
      (assert-same-copy 'lla-integer lla::+single+)
      (assert-same-copy 'lla-single lla::+single+)
      (assert-same-copy 'lla-double lla::+single+)
      (assert-copy-error 'lla-complex-single lla::+single+)
      ;; double float
      (assert-same-copy t lla::+double+)
      (assert-same-copy 'lla-integer lla::+double+)
      (assert-same-copy 'lla-single lla::+double+)
      (assert-same-copy 'lla-double lla::+double+)
      (assert-copy-error 'lla-complex-single lla::+double+)
      ;; complex single
      (assert-same-copy t lla::+complex-single+)
      (assert-same-copy 'lla-integer lla::+complex-single+)
      (assert-same-copy 'lla-single lla::+complex-single+)
      (assert-same-copy 'lla-double lla::+complex-single+)
      (assert-same-copy 'lla-complex-single lla::+complex-single+)
      (assert-same-copy 'lla-complex-double lla::+complex-single+)
      ;; complex double
      (assert-same-copy t lla::+complex-double+)
      (assert-same-copy 'lla-integer lla::+complex-double+)
      (assert-same-copy 'lla-single lla::+complex-double+)
      (assert-same-copy 'lla-double lla::+complex-double+)
      (assert-same-copy 'lla-complex-single lla::+complex-double+)
      (assert-same-copy 'lla-complex-double lla::+complex-double+))))

(deftest transpose-to-memory (pinned-array-suite)
  (with-foreign-temporary-buffer (pointer 100)
    (let+ (((&flet check-transpose (matrix internal-type)
              (lla::transpose-matrix-to-memory matrix pointer internal-type)
              (assert-true (check-memory-contents2 pointer internal-type
                                                   (aops:transpose matrix)))))
           ((&flet assert-same-transpose (element-type internal-type)
              (check-transpose (random-array element-type 5 6)
                               internal-type)))
           ((&flet assert-transpose-error (element-type internal-type)
              (assert-condition error
                  (lla::transpose-matrix-to-memory
                   (random-array element-type 7 5) pointer internal-type)))))
      ;; single float
      (assert-same-transpose t lla::+single+)
      (assert-same-transpose 'lla-integer lla::+single+)
      (assert-same-transpose 'lla-single lla::+single+)
      (assert-same-transpose 'lla-double lla::+single+)
      (assert-transpose-error 'lla-complex-single lla::+single+)
      ;; ;; double float
      (assert-same-transpose t lla::+double+)
      (assert-same-transpose 'lla-integer lla::+double+)
      (assert-same-transpose 'lla-single lla::+double+)
      (assert-same-transpose 'lla-double lla::+double+)
      (assert-transpose-error 'lla-complex-single lla::+double+)
      ;; complex single
      (assert-same-transpose t lla::+complex-single+)
      (assert-same-transpose 'lla-integer lla::+complex-single+)
      (assert-same-transpose 'lla-single lla::+complex-single+)
      (assert-same-transpose 'lla-double lla::+complex-single+)
      (assert-same-transpose 'lla-complex-single lla::+complex-single+)
      (assert-same-transpose 'lla-complex-double lla::+complex-single+)
      ;; complex double
      (assert-same-transpose t lla::+complex-double+)
      (assert-same-transpose 'lla-integer lla::+complex-double+)
      (assert-same-transpose 'lla-single lla::+complex-double+)
      (assert-same-transpose 'lla-double lla::+complex-double+)
      (assert-same-transpose 'lla-complex-single lla::+complex-double+)
      (assert-same-transpose 'lla-complex-double lla::+complex-double+))))

(deftest copy-from-memory (pinned-array-suite)
  (with-foreign-temporary-buffer (pointer 100)
    (let+ (((&flet check-copy-to-and-from (element-type internal-type)
              (let ((array (random-array element-type 5 8)))
                (lla::copy-array-to-memory array pointer internal-type)
                (assert-equalp array
                    (lla::create-array-from-memory
                     pointer internal-type (array-dimensions array)))))))
      (check-copy-to-and-from 'lla-single lla::+single+)
      (check-copy-to-and-from 'lla-double lla::+double+)
      (check-copy-to-and-from 'lla-complex-single lla::+complex-single+)
      (check-copy-to-and-from 'lla-complex-double lla::+complex-double+)
      (check-copy-to-and-from 'lla-integer lla::+integer+))))

(deftest transpose-from-memory (pinned-array-suite)
  (with-foreign-temporary-buffer (pointer 100)
    (let+ (((&flet check-transpose-to-and-from (element-type internal-type)
              (let ((matrix (random-array element-type 5 8)))
                (lla::transpose-matrix-to-memory matrix pointer internal-type)
                (assert-equalp matrix
                    (lla::create-transposed-matrix-from-memory
                     pointer internal-type (array-dimensions matrix)))))))
      (check-transpose-to-and-from 'lla-single lla::+single+)
      (check-transpose-to-and-from 'lla-double lla::+double+)
      (check-transpose-to-and-from 'lla-complex-single lla::+complex-single+)
      (check-transpose-to-and-from 'lla-complex-double lla::+complex-double+)
      (check-transpose-to-and-from 'lla-integer lla::+integer+))))

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

(deftest (pinned-array-suite)
  pinning-input
  (iter
    (for (source destination) :in (coercible-pairs-list))
    (ensure
     (test-pinning-readonly source destination))))

;; pinning, copy

(deftest (pinned-array-suite)
  pinning-copy
  (iter
    (for (source destination) :in (coercible-pairs-list))
    (ensure
     (test-pinning-copy source destination))))

;; pinning, output

(deftest (pinned-array-suite)
  pinning-output
  (iter
    (for (source destination) :in (coercible-pairs-list))
    (ensure
     (test-pinning-output source destination))))

;; vector-output

(deftest (pinned-array-suite)
  vector-output
  (ensure (test-vector-output lla::+integer+))
  (ensure (test-vector-output lla::+single+))
  (ensure (test-vector-output lla::+double+))
  (ensure (test-vector-output lla::+complex-single+))
  (ensure (test-vector-output lla::+complex-double+)))
