(in-package :lla)

;;;; Interface
;;;
;;; The following macros provide the abstract interface used by the
;;; rest of LLA. When implementation-specific speedups (eg shared
;;; arrays) are not available, they fall back to copying.

(defmacro with-array-input (((pointer &optional (copied? (gensym "COPIED?")))
                             array internal-type transpose? force-copy?)
                            &body body)
  "Ensure that ARRAY is mapped to a corresponding memory area for the duration of BODY (see below for details of semantics).  POINTER is bound to the start of the memory address.  The representation of values in the memory is determined by INTERNAL-TYPE.  When TRANSPOSE?, transpose the array before BODY (only works for matrices, otherwise signal an error).

If FORCE-COPY? is false, POINTER may or may not point to the memory backing ARRAY.

If FORCE-COPY? is true, POINTER always points to a copy of the contents of ARRAY.

COPIED? is bound to indicate whether POINTER points to a copy or the actual array contents.

The value of the expression is always the value of BODY."
  (alexandria:with-gensyms (backing-array offset)
    (once-only (array internal-type transpose?)
      `(flet ((body (,pointer ,copied?)
                (declare (ignorable ,copied?))
                ,@body))
         (declare (dynamic-extent #'body)
                  #+sbcl (sb-ext:muffle-conditions sb-ext:code-deletion-note))
         (multiple-value-bind (,backing-array ,offset) (backing-array ,array)
           (if (and (not ,transpose?)
                    (not ,force-copy?)
                    (shareable-array? ,backing-array ,internal-type))
               (with-pinned-array (,pointer ,backing-array)
                 (cffi:incf-pointer ,pointer
                     (* ,offset (foreign-size ,internal-type)))
                 (body ,pointer nil))
               (with-work-area (,pointer ,internal-type
                                         (array-total-size ,array))
                 (if ,transpose?
                     (transpose-matrix-to-memory ,array ,pointer ,internal-type)
                     (copy-array-to-memory ,array ,pointer ,internal-type))
                 (body ,pointer t))))))))

(defmacro with-array-output (((pointer &optional (copied? (gensym "COPIED?")))
                              array internal-type transpose?)
                             &body body)
  "Ensure that ARRAY is mapped to a corresponding memory area for the duration of BODY (see below for details of semantics).  POINTER is bound to the start of the memory address.  The representation of values in the memory is determined by INTERNAL-TYPE.  When TRANSPOSE?, transpose the array after BODY (only works for matrices, otherwise signal an error).

COPIED? is bound to indicate whether POINTER points to a copy or the actual array contents.

The value of the expression is always the value of BODY."
  (alexandria:with-gensyms (backing-array offset)
    (once-only (array internal-type)
      `(flet ((body (,pointer ,copied?)
                (declare (ignorable ,copied?))
                ,@body))
         (declare (dynamic-extent #'body)
                  #+sbcl (sb-ext:muffle-conditions sb-ext:code-deletion-note))
         (multiple-value-bind (,backing-array ,offset) (backing-array ,array)
           (if (and (not ,transpose?)
                    (shareable-array? ,backing-array ,internal-type))
               (with-pinned-array (,pointer ,backing-array)
                 (cffi:incf-pointer ,pointer
                     (* ,offset (foreign-size ,internal-type)))
                 (body ,pointer nil))
               (with-work-area (,pointer ,internal-type
                                         (array-total-size ,array))
                 (multiple-value-prog1 (body ,pointer t)
                   (if ,transpose?
                       (transpose-matrix-from-memory ,array ,pointer
                                                     ,internal-type)
                       (copy-array-from-memory ,array ,pointer
                                               ,internal-type))))))))))

(defmacro with-array-input-output (((pointer
                                     &optional (copied? (gensym "COPIED?")))
                                    input-array input-internal-type
                                    input-transpose? input-force-copy?
                                    output-array output-internal-type
                                    output-transpose?)
                                   &body body)
  "Similar to WITH-ARRAY-INPUT, it also ensures that OUTPUT-ARRAY contains the contents the memory area pointed to by POINTER when BODY is finished. If  OUTPUT-ARRAY is NIL then it is equivalent to WITH-ARRAY-INPUT."
  (once-only (input-array output-array)
    `(with-array-input ((,pointer ,copied?) ,input-array
                        ,input-internal-type ,input-transpose?
                        ,input-force-copy?)
       (multiple-value-prog1 (progn ,@body)
         (locally
             (declare
              #+sbcl (sb-ext:muffle-conditions sb-ext:code-deletion-note))
           (cond ((null ,output-array)
                  ;; The output is unused.
                  nil)
                 ((and (eq ,input-array ,output-array)
                       (not ,copied?))
                  (when ,output-transpose?
                    (transpose-matrix-from-memory ,output-array ,pointer
                                                  ,output-internal-type)))
                 (t
                  (if ,output-transpose?
                      (transpose-matrix-from-memory ,output-array ,pointer
                                                    ,output-internal-type)
                      (copy-array-from-memory ,output-array ,pointer
                                              ,output-internal-type)))))))))

(defmacro with-work-area ((pointer internal-type size) &body body)
  "Allocate a work area of SIZE and INTERNAL-TYPE, and bind the POINTER to its start during BODY."
  (check-type pointer symbol)
  `(with-foreign-pointer (,pointer (* (foreign-size ,internal-type) ,size))
     ,@body))

;;;; Implementation
;;;
;;; We copy the array to and from the memory area if necessary, either because it must be transposed, copying is requested explicitly, true pinning is not implemented or LLA:CFFI-PINNING is enabled in the configuration.

(defun backing-array (array)
  "Return the array in which the contents of ARRAY are stored. For simple arrays, this is always the array itself.  The second value is the displacement."
  #+sbcl
  (sb-c::with-array-data ((v array) (start) (end))
    (declare (ignore end))
    (values v start))
  #-sbcl
  (values array 0))

#+(and sbcl (not lla::cffi-pinning))
(defun shareable-array? (array internal-type)
  (and (equal (upgraded-array-element-type (lisp-type internal-type))
              (array-element-type array))
       (typep array 'simple-array)))

#-(and sbcl (not lla::cffi-pinning))
(defun shareable-array? (array internal-type)
  (declare (ignore array internal-type))
  nil)

#+sbcl
(defmacro with-pinned-array ((pointer array) &body body)
  (once-only (array)
    `(sb-sys:with-pinned-objects (,array)
       (let ((,pointer (sb-sys:vector-sap
                        (sb-ext:array-storage-vector ,array))))
         ,@body))))

#-sbcl
(defmacro with-pinned-array ((pointer array) &body body)
  (declare (ignore pointer array body))
  `(error "Pinning is not implemented on this platform."))
