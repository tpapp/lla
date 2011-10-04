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

(defun check-memory-contents2 (pointer lla-type values)
  "Like CHECK-MEMORY-CONTENTS, with autoconversion of values and calculation
of LLA-TYPE."
  (check-memory-contents pointer
                         (ecase lla-type
                           ((:single :complex-single) :float)
                           ((:double :complex-double) :double)
                           (:integer #-lla::int64 :int32
                                     #+lla::int64 :int64))
                         (let ((vector (aetypecase values
                                         (vector it)
                                         (array (flatten-array it))
                                         (t (coerce it 'vector)))))
                           (case lla-type
                             ((:complex-single :complex-double)
                              (flatten-complex-vector vector))
                             (t vector)))))

(addtest (pinned-array-tests)
  copy-to-memory
  (with-foreign-temporary-buffer (pointer 100)
    (let+ (((&flet check-copy (array lla-type)
              (lla::copy-array-to-memory array pointer lla-type)
              (check-memory-contents2 pointer lla-type array)))
           ((&macrolet ensure-copy (element-type lla-type)
              `(ensure (check-copy (random-array ,element-type 5 7)
                                   ,lla-type))))
           ((&macrolet ensure-copy-error (element-type lla-type)
              `(ensure-error
                 (lla::copy-array-to-memory (random-array ,element-type 19)
                                            pointer ,lla-type)))))
      ;; single float
      (ensure-copy t :single)
      (ensure-copy 'lla-integer :single)
      (ensure-copy 'lla-single :single)
      (ensure-copy 'lla-double :single)
      (ensure-copy-error 'lla-complex-single :single)
      ;; double float
      (ensure-copy t :double)
      (ensure-copy 'lla-integer :double)
      (ensure-copy 'lla-single :double)
      (ensure-copy 'lla-double :double)
      (ensure-copy-error 'lla-complex-single :double)
      ;; complex single
      (ensure-copy t :complex-single)
      (ensure-copy 'lla-integer :complex-single)
      (ensure-copy 'lla-single :complex-single)
      (ensure-copy 'lla-double :complex-single)
      (ensure-copy 'lla-complex-single :complex-single)
      (ensure-copy 'lla-complex-double :complex-single)
      ;; complex double
      (ensure-copy t :complex-double)
      (ensure-copy 'lla-integer :complex-double)
      (ensure-copy 'lla-single :complex-double)
      (ensure-copy 'lla-double :complex-double)
      (ensure-copy 'lla-complex-single :complex-double)
      (ensure-copy 'lla-complex-double :complex-double))))

(addtest (pinned-array-tests)
  transpose-to-memory
  (with-foreign-temporary-buffer (pointer 100)
    (let+ (((&flet check-transpose (matrix lla-type)
              (lla::transpose-matrix-to-memory matrix pointer lla-type)
              (check-memory-contents2 pointer lla-type (transpose matrix))))
           ((&macrolet ensure-transpose (element-type lla-type)
              `(ensure (check-transpose (random-array ,element-type 5 6)
                                        ,lla-type))))
           ((&macrolet ensure-transpose-error (element-type lla-type)
              `(ensure-error 
                 (lla::transpose-matrix-to-memory 
                  (random-array ,element-type 7 5) pointer ,lla-type)))))
      ;; single float
      (ensure-transpose t :single)
      (ensure-transpose 'lla-integer :single)
      (ensure-transpose 'lla-single :single)
      (ensure-transpose 'lla-double :single)
      (ensure-transpose-error 'lla-complex-single :single)
      ;; ;; double float
      (ensure-transpose t :double)
      (ensure-transpose 'lla-integer :double)
      (ensure-transpose 'lla-single :double)
      (ensure-transpose 'lla-double :double)
      (ensure-transpose-error 'lla-complex-single :double)
      ;; complex single
      (ensure-transpose t :complex-single)
      (ensure-transpose 'lla-integer :complex-single)
      (ensure-transpose 'lla-single :complex-single)
      (ensure-transpose 'lla-double :complex-single)
      (ensure-transpose 'lla-complex-single :complex-single)
      (ensure-transpose 'lla-complex-double :complex-single)
      ;; complex double
      (ensure-transpose t :complex-double)
      (ensure-transpose 'lla-integer :complex-double)
      (ensure-transpose 'lla-single :complex-double)
      (ensure-transpose 'lla-double :complex-double)
      (ensure-transpose 'lla-complex-single :complex-double)
      (ensure-transpose 'lla-complex-double :complex-double))))

(addtest (pinned-array-tests)
  copy-from-memory
  (with-foreign-temporary-buffer (pointer 100)
    (let+ (((&flet check-copy-to-and-from (element-type lla-type)
              (let ((array (random-array element-type 5 8)))
                (lla::copy-array-to-memory array pointer lla-type)
                (equalp array (lla::create-array-from-memory
                               pointer lla-type (array-dimensions array))))))
           ((&macrolet ensure-to-and-from (element-type lla-type)
              `(ensure (check-copy-to-and-from ,element-type ,lla-type)))))
      (ensure-to-and-from 'lla-single :single)
      (ensure-to-and-from 'lla-double :double)
      (ensure-to-and-from 'lla-complex-single :complex-single)
      (ensure-to-and-from 'lla-complex-double :complex-double)
      (ensure-to-and-from 'lla-integer :integer))))

(addtest (pinned-array-tests)
  transpose-from-memory
  (with-foreign-temporary-buffer (pointer 100)
    (let+ (((&flet check-transpose-to-and-from (element-type lla-type)
              (let ((matrix (random-array element-type 5 8)))
                (lla::transpose-matrix-to-memory matrix pointer lla-type)
                (equalp matrix (lla::create-transposed-matrix-from-memory
                                pointer lla-type (array-dimensions matrix))))))
           ((&macrolet ensure-to-and-from (element-type lla-type)
              `(ensure (check-transpose-to-and-from ,element-type ,lla-type)))))
      (ensure-to-and-from 'lla-single :single)
      (ensure-to-and-from 'lla-double :double)
      (ensure-to-and-from 'lla-complex-single :complex-single)
      (ensure-to-and-from 'lla-complex-double :complex-double)
      (ensure-to-and-from 'lla-integer :integer))))

;;;; Pinning tests
;;; 
;;; Test the pinning framework, regardless of its implementation (copying,
;;; sharing, etc).

(defun coercible? (lla-source-type lla-target-type)
  "Permitted coercions for LLA types.  It is guaranteed that CL:COERCE can perform
these coercions on the corresponding lisp types. Basic summary: (1) integers can be
coerced to anything, (2) single<->double precision coercions are possible both
ways, (3) real floats can be upgraded to complex."
  (cond
    ;; always valid
    ((eq lla-source-type :integer) t)
    ;; nothing else can be converted to integer
    ((eq lla-target-type :integer) nil)
    ;; no complex->real
    ((and (lla-complex? lla-source-type)
          (not (lla-complex? lla-target-type)))
     nil)
    ;; all the rest should be possible
    (t t)))

(defun coercible-pairs-list ()
  ;;   "Generate the list of all LLA (source target) pairs for which
  ;; coercible? holds.  For internal use only, NOT EXPORTED."
  (let (coercible)
    (dolist (source (lla-types))
      (dolist (target (lla-types))
        (when (coercible? source target)
          (push (list source target) coercible))))
    ;; reverse only for cosmetic purposes (debugging)
    (nreverse coercible)))

(defun test-pinning-readonly (source-type destination-type &key
                              (report? t) (length 50))
  (let ((vector (make-random-array length source-type 100)))
    (lla::with-pinned-array (pointer vector destination-type nil nil nil nil)
      (array-at-pointer= vector pointer destination-type :report? report?))))

(defun test-pinning-copy (source-type destination-type &key
                          (report? t) (length 50))
  (let* ((vector (make-random-array length source-type 100))
         (copy (copy-seq vector))
         (copy-inc (map 'vector #'1+ copy)))
    (lla::with-pinned-array (pointer vector destination-type nil
                                     :copy nil nil)
      ;; check equality
      (unless (array-at-pointer= copy pointer destination-type
                                 :report? report?)
        (return-from test-pinning-copy nil))
      ;; increment memory area, check equality, and that original is intact
      (dotimes (i length)
        (incf (lla::mem-aref* pointer destination-type i) 1))
      (and (array-at-pointer= copy-inc pointer destination-type
                              :report? report?)
           (equalp vector copy)))))

(defun test-pinning-output (source-type destination-type &key
                            (report? t) (length 50))
  (let* ((vector (make-random-array length source-type 100))
         vector-inc
         (copy (copy-seq vector))
         (copy-inc (map 'vector #'1+ copy)))
    (lla::with-pinned-array (pointer vector destination-type nil
                                    vector-inc nil nil)
      ;; check equality
      (unless (array-at-pointer= copy pointer destination-type :report? report?)
        (return-from test-pinning-output nil))
      ;; increment memory area, check equality, and that original is intact
      (dotimes (i length)
        (incf (lla::mem-aref* pointer destination-type i) 1))
      (unless (and (array-at-pointer= copy-inc pointer destination-type
                                      :report? report?)
                   (equalp vector copy))
        (return-from test-pinning-output nil)))
    ;; check output
    (every #'= vector-inc copy-inc)))

(defun test-vector-output (lla-type &key (length 50))
  (let ((vector (make-random-array length lla-type 100))
        output)
    (lla::with-array-output (pointer output lla-type length nil)
      (dotimes (index length)
        (setf (lla::mem-aref* pointer lla-type index) (aref vector index))))
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
  (ensure (test-vector-output :integer))
  (ensure (test-vector-output :single))
  (ensure (test-vector-output :double))
  (ensure (test-vector-output :complex-single))
  (ensure (test-vector-output :complex-double)))


