;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla-tests)

(deftestsuite pinned-array-tests (lla-tests)
  ())

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


