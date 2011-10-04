(in-package :lla)

(defmacro with-pinned-array ((pointer array lla-type transpose? 
                              output output-dimensions output-transpose?)
                             &body body)
  "Ensure that ARRAY is mapped to a corresponding memory area for the duration
  of BODY (see below for details of semantics).  POINTER is bound to the start
  of the memory address.  The representation of values in the memory is
  determined by LLA-TYPE.  When TRANSPOSE?, transpose the array before
  BODY (only works for matrices, otherwise signal an error).  TRANSPOSE? has
  to be fixed at compile-time.

  The semantics of memory access is determined by OUTPUT.

  1. If OUTPUT is NIL, the memory area is NOT expected to be modified by BODY.

  2. If OUTPUT is :COPY, the memory area may be modified by BODY.

  3. If OUTPUT is anything else, the memory area may be modifed by BODY, and
  the result will be available after that while BODY is executed, assigned to
  the place referred to by OUTPUT, with the given OUTPUT-DIMENSIONS (which
  default to the dimensions of ARRAY when NIL).  When OUTPUT-TRANSPOSE?, a
  transposed result is saved.

  The value of the expression is always the value of body."
  `(#+lla::cffi-pinning with-pinned-array-copy%
    #+(and sbcl (not lla::cffi-pinning)) with-pinned-array-sbcl%
    (,pointer ,array ,lla-type ,transpose? 
              ,output ,output-dimensions ,output-transpose?)
    ,@body))

(defmacro with-array-output ((pointer output lla-type dimensions transpose?)
                             &body body)
  "Allocate a memory area of the given LLA-TYPE and DIMENSIONS, bind its
  address to POINTER, and copy the contents to OUTPUT after BODY.  When the
  implied total size is 0, the implementation is allowed to assign the null
  pointer to POINTER.  When TRANSPOSE?, matrices are transposed."
  `(#+lla::cffi-pinning with-array-output-copy%
    #+(and sbcl (not lla::cffi-pinning)) with-array-output-sbcl%
    (,pointer ,output ,lla-type ,dimensions ,transpose?)
    ,@body))


;;; Pinning implemented by copying.
;;;
;;; We effectively simulate pinning by copying the array to and from the
;;; memory area.

(defmacro with-copied-elements% ((pointer array lla-type transpose?)
                                 &body body)
  "Allocate memory and copy elements of ARRAY to POINTER for the scope of
BODY.  When TRANSPOSE?, matrices are transposed."
  (check-type pointer symbol)
  (check-type transpose? boolean)
  (once-only (array lla-type)
    (with-unique-names (length size)
      `(let ((,length (array-total-size ,array))
             (,size (foreign-size* ,lla-type)))
         (with-foreign-pointer (,pointer (* ,size ,length))
           ,(if transpose?
                `(transpose-matrix-to-memory ,array ,pointer ,lla-type)
                `(copy-array-to-memory ,array ,pointer ,lla-type))
           ,@body)))))

(defmacro with-pinned-array-copy% ((pointer array lla-type transpose? output
                                    output-dimensions output-transpose?)
                                   &body body)
  ;; note: OUTPUT-DIMENSIONS used only once, and thus not wrapped in ONCE-ONLY
  (check-type pointer symbol)
  (check-type output-transpose? boolean)
  (once-only (array lla-type)
    `(with-copied-elements% (,pointer ,array ,lla-type ,transpose?)
       (multiple-value-prog1
           (progn ,@body)
         ,@(unless (or (not output) (eq output :copy))
             `((setf ,output
                     ,(let ((output-dimensions
                              (aif output-dimensions
                                   it
                                   `(array-dimensions ,array))))
                        (if output-transpose?
                            `(create-transposed-matrix-from-memory
                              ,pointer ,lla-type ,output-dimensions)
                            `(create-array-from-memory
                              ,pointer ,lla-type ,output-dimensions))))))))))

(defmacro with-array-output-copy% ((pointer output lla-type dimensions
                                    transpose?)
                                   &body body)
  (check-type pointer symbol)
  (once-only (lla-type dimensions)
    `(with-foreign-pointer (,pointer
                            (* (foreign-size* ,lla-type)
                               (reduce #'* (ensure-list ,dimensions))))
       (multiple-value-prog1
           (progn ,@body)
         (setf ,output
               ,(if transpose?
                    `(create-transposed-matrix-from-memory ,pointer ,lla-type
                                                           ,dimensions)
                    `(create-array-from-memory ,pointer ,lla-type
                                               ,dimensions)))))))

(defmacro with-work-area ((pointer lla-type size) &body body)
  (check-type pointer symbol)
  `(with-foreign-pointer (,pointer (* (foreign-size* ,lla-type) ,size))
     ,@body))

;;; SBCL specific code

;; #+sbcl
;; (defun ensure-sharable-array (array lla-type copy?)
;;   ""
;;   (if (and (equal (lla-to-lisp-type lla-type) (array-element-type array))
;;            (typep array 'simple-array)
;;            (not copy?))
;;       array
;;       (convert-array* array lla-type)))

;; #+sbcl
;; (defmacro with-pinned-array-sbcl% (pointer array lla-type transpose? output
;;                                    output-dimensions output-transpose?
;;                                    &body body)
;;   (let ((array-var (gensym))
;;         (copy? (when output t)))
;;     `(let ((,array-var (ensure-sharable-array ,array ,lla-type ,copy?)))
;;        (sb-sys:with-pinned-objects (,array-var)
;;          (let ((,pointer (sb-sys:vector-sap
;;                           (sb-ext:array-storage-vector ,array-var))))
;;            ,@body
;;            ,@(when (and output (not (eq output :copy)))
;;                `((setf ,output ,array-var))))))))





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
