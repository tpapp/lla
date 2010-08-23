;;; -*- Mode:Lisp; Syntax:ANSI-Common-Lisp; Coding:utf-8 -*-

(in-package #:lla-unit-tests)

(deftestsuite pinned-vector-tests (lla-unit-tests)
  ())

(defun coercible-p (lla-source-type lla-target-type)
  "Permitted coercions for LLA types.  It is guaranteed that
CL:COERCE can perform these coercions on the corresponding lisp
types. Basic summary: (1) integers can be coerced to anything, \(2)
single<->double precision coercions are possible both ways, (3) real
floats can be upgraded to complex."
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
  ;; coercible-p holds.  For internal use only, NOT EXPORTED."
  (let ((lla-types '(:integer :single :double :complex-single :complex-double))
        coercible)
    (dolist (source lla-types)
      (dolist (target lla-types)
        (when (coercible-p source target)
          (push (list source target) coercible))))
    ;; reverse only for cosmetic purposes
    (nreverse coercible)))

(defun test-pinning-readonly (source-type destination-type &key
                              (report-p t) (length 50))
  (let ((vector (make-random-vector length source-type 100)))
    (with-pinned-vector (vector pointer destination-type)
      (vector-at-pointer= vector pointer destination-type :report-p
                          report-p))))

(defun test-pinning-copy (source-type destination-type &key
                          (report-p t) (length 50))
  (let* ((vector (make-random-vector length source-type 100))
         (copy (copy-seq vector))
         (copy-inc (map 'vector #'1+ copy)))
    (with-pinned-vector (vector pointer destination-type :copy)
      ;; check equality
      (unless (vector-at-pointer= copy pointer destination-type :report-p report-p)
        (return-from test-pinning-copy nil))
      ;; increment memory area, check equality, and that original is intact
      (dotimes (i length)
        (incf (mem-aref* pointer destination-type i) 1))
      (and (vector-at-pointer= copy-inc pointer destination-type
                               :report-p report-p)
           (equalp vector copy)))))

(defun test-pinning-output (source-type destination-type &key
                            (report-p t) (length 50))
  (let* ((vector (make-random-vector length source-type 100))
         vector-inc
         (copy (copy-seq vector))
         (copy-inc (map 'vector #'1+ copy)))
    (with-pinned-vector (vector pointer destination-type vector-inc)
      ;; check equality
      (unless (vector-at-pointer= copy pointer destination-type :report-p report-p)
        (return-from test-pinning-output nil))
      ;; increment memory area, check equality, and that original is intact
      (dotimes (i length)
        (incf (mem-aref* pointer destination-type i) 1))
      (unless (and (vector-at-pointer= copy-inc pointer destination-type
                                       :report-p report-p)
                   (equalp vector copy))
        (return-from test-pinning-output nil)))
    ;; check output
    (every #'= vector-inc copy-inc)))

(defun test-vector-output (lla-type &key (length 50))
  (let ((vector (make-random-vector length lla-type 100))
        output)
    (lla::with-vector-output (output pointer lla-type length)
      (dotimes (index length)
        (setf (mem-aref* pointer lla-type index) (aref vector index))))
    (equalp vector output)))


;; pinning, input only

(addtest (pinned-vector-tests)
  pinning-input
  (iter 
    (for (source destination) :in (coercible-pairs-list))
    (ensure
     (test-pinning-readonly source destination))))

;; pinning, copy

(addtest (pinned-vector-tests)
  pinning-copy
  (iter 
    (for (source destination) :in (coercible-pairs-list))
    (ensure
     (test-pinning-copy source destination))))

;; pinning, output

(addtest (pinned-vector-tests)
  pinning-output
  (iter 
    (for (source destination) :in (coercible-pairs-list))
    (ensure
     (test-pinning-output source destination))))

;; vector-output

(addtest (pinned-vector-tests)
  vector-output
  (ensure (test-vector-output :integer))
  (ensure (test-vector-output :single))
  (ensure (test-vector-output :double))
  (ensure (test-vector-output :complex-single))
  (ensure (test-vector-output :complex-double)))


