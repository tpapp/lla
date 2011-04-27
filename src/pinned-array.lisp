(in-package :lla)

;;; The WITH-PINNED-ARRAY macro maps an array to a memory area of values with
;;; specified LLA-TYPE, either for input, input-output or output, converting the
;;; input if required.

;;; If :PINNED-ARRAY-COPY is in *features*, the macros will use copying to
;;; emulate pinned vectors regardless of the implementation.

(eval-when (:compile-toplevel :load-toplevel :execute)
  ;; #-sbcl 
  (pushnew :pinned-array-copy *features*))

;;; All-in-one wrapper macro for numeric vectors.  Call one of the
;;; implementation-specific macros below, should be extended accordingly.

(defmacro with-pinned-array ((array pointer lla-type &rest rest &key output
                                    output-dimensions)
                             &body body)
  "Ensure that ARRAY is mapped to a corresponding memory area for the duration
  of BODY (see below for details of semantics).  POINTER is bound to the start
  of the memory address.  The representation of values in the memory is
  determined by LLA-TYPE.

  The semantics is determined by OUTPUT.

  1. If OUTPUT is NIL (or not given), the memory area is not expected to be
  modified by BODY.

  2. If OUTPUT is :COPY, the memory area may be modified by BODY.

  3. If OUTPUT is anything else, the memory area may be modifed by BODY, and the
  result will be available after that while BODY is executed, assigned to the place
  referred to by output, with the given OUTPUT-DIMENSIONS (which default to the
  dimensions of the array).

  The value of the expression is always the value of body."
  (declare (ignorable output output-dimensions))
  `(#+pinned-array-copy with-pinned-array-copy%
    #+(and sbcl (not pinned-array-copy)) with-pinned-array-sbcl%
    (,array ,pointer ,lla-type ,@rest)
    ,@body))

(define-with-multiple-bindings with-pinned-array)

(defmacro with-array-output ((output pointer lla-type dimensions)
                             &body body)
  "Allocate a memory area of the given LLA-TYPE and DIMENSIONS, bind its
  address to POINTER, and copy the contents to OUTPUT after BODY.  When the
  implied total size is 0, the implementation is allowed to assign the null
  pointer to POINTER."
  `(#+pinned-array-copy with-array-output-copy%
    #+(and sbcl (not pinned-array-copy)) with-array-output-sbcl%
    (,output ,pointer ,lla-type ,dimensions)
    ,@body))

(define-with-multiple-bindings with-array-output)


;;; Pinning implemented by copying.
;;;
;;; We effectively simulate pinning by copying the array to and from the memory
;;; area.

(defmacro with-copied-elements% ((array pointer lla-type) &body body)
  (check-type pointer symbol)
  (once-only (array lla-type)
    (with-unique-names (length size index)
      `(let ((,length (array-total-size ,array))
             (,size (foreign-size* ,lla-type)))
         (with-foreign-pointer (,pointer (* ,size ,length))
           (dotimes (,index ,length)
             (setf (mem-aref* ,pointer ,lla-type ,index)
                   (coerce* (row-major-aref ,array ,index) ,lla-type)))
           ,@body)))))

(defun copy-array-from-memory% (pointer lla-type dimensions)
  ;; !! this could be speeded up by conditioning on lla-type
  (let ((array (make-array* dimensions lla-type)))
    (dotimes (index (array-total-size array))
      (setf (row-major-aref array index)
            (mem-aref* pointer lla-type index)))
    array))

(defmacro with-pinned-array-copy% ((array pointer lla-type &key
                                          output output-dimensions)
                                   &body body)
  ;; note: OUTPUT-DIMENSIONS used only once, and thus not wrapped in ONCE-ONLY
  (check-type pointer symbol)
  (once-only (array lla-type)
    `(with-copied-elements% (,array ,pointer ,lla-type)
       (multiple-value-prog1
           (progn ,@body)
         ,@(unless (or (not output) (eq output :copy))
             `((setf ,output 
                     (copy-array-from-memory% ,pointer ,lla-type
                                              ,(aif output-dimensions
                                                    output-dimensions
                                                    `(array-dimensions ,array))))))))))

(defmacro with-array-output-copy% ((output pointer lla-type dimensions)
                                   &body body)
  (check-type pointer symbol)
  (once-only (lla-type dimensions)
    `(with-foreign-pointer (,pointer (* (foreign-size* ,lla-type)
                                        (reduce #'* (ensure-list ,dimensions))))
       (multiple-value-prog1
           (progn ,@body)
         (setf ,output
               (copy-array-from-memory% ,pointer ,lla-type ,dimensions))))))

;;; SBCL specific code

;; #+sbcl
;; (defmacro pinned-vector-wrapper-sbcl% ((vector pointer) &body body)
;;   "Pin the vector and bind pointer to its data during body.  This is a
;; utility function for implementations with pinning, and should not be
;; used directly in other files, ie it is not part of the interface."
;;   (check-type pointer symbol)
;;   (once-only (vector)
;;     `(sb-sys:with-pinned-objects (,vector)
;;        (let ((,pointer (sb-sys:vector-sap ,vector)))
;; 	 ,@body))))

;; #+sbcl
;; (defmacro with-pinned-vector-sbcl% ((vector pointer foreign-type &optional output)
;;                                     &body body)
;;   (once-only (foreign-type vector)
;;     (if output
;;         (let ((vector-copy (gensym* pointer '#:-copy)))
;;           `(let ((,vector-copy (copy-vector ,vector ,foreign-type)))
;;              (check-type ,foreign-type foreign-type)
;;              (multiple-value-prog1 
;;                  (pinned-vector-wrapper-sbcl% (,vector-copy ,pointer)
;;                    ,@body)
;;                ,@(unless (eq output :copy)
;;                    `((setf ,output ,vector-copy))))))
;;         `(pinned-vector-wrapper-sbcl% 
;;              ((if (and (simple-array1? ,vector)
;;                        (eq (array-foreign-type ,vector)
;;                            ,foreign-type))
;;                   ,vector
;;                   (copy-vector ,vector ,foreign-type))
;;               ,pointer)
;;            ,@body))))

;; #+sbcl
;; (defmacro with-vector-output-sbcl% ((output pointer foreign-type length)
;;                                     &body body)
;;   (check-type pointer symbol)
;;   (once-only (foreign-type length)
;;     (let ((vector (gensym* pointer)))
;;       `(let ((,vector (lla-array ,length ,foreign-type)))
;;          (check-type ,foreign-type foreign-type)
;;          (multiple-value-prog1
;;              (pinned-vector-wrapper-sbcl% (,vector ,pointer)
;;                ,@body)
;;            (setf ,output ,vector))))))
