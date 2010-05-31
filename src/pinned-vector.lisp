(in-package :lla)

;;; Pinned vector macros map a numeric vector to a memory area of
;;; value wich specified LLA-TYPE, either for input, input-output or
;;; output, converting the input if required.

;;; If :lla-pinned-vector-cffi is in *features*, LLA will use CFFI to emulate
;;; pinned vectors.

(eval-when (:compile-toplevel :load-toplevel :execute)
  #-sbcl (pushnew :lla-pinned-vector-cffi *features*))

;;; All-in-one wrapper macro for numeric vectors.  Call one of the
;;; implementation-specific macros below, should be extended accordingly.

(defmacro with-pinned-vector ((vector pointer lla-type &optional output)
                              &body body)
  "Ensure that vector is mapped to a corresponding memory area for the
  duration of BODY (see below for details of semantics).  POINTER is
  bound to the start of the memory address.  The representation of
  values in the memory is determined by LLA-TYPE.

  The semantics is determined by OUTPUT.

  1. If OUTPUT is NIL (or not given), the memory area may not be
  modified by BODY.

  2. If OUTPUT is :COPY, the memory area may be modified by BODY.

  3. If OUTPUT is anything else, the memory area may be modifed by
  BODY, and the result will be available after that while BODY is
  executed, assigned to the place referred to by output.

  The value of the expression is always the value of body."
  `(#+lla-pinned-vector-cffi with-pinned-vector-general% 
    #+sbcl with-pinned-vector-sbcl%
    (,vector ,pointer ,lla-type ,output)
    ,@body))

(define-with-multiple-bindings with-pinned-vector)

(defmacro with-vector-output ((output pointer lla-type length)
                              &body body)
  "Allocate a memory area of the given LLA-TYPE and LENGTH, bind its address to
  POINTER, and copy the contents to OUTPUT after BODY.  When (ZEROP LENGTH), the
  implementation is allowed to assign the null pointer to pointer."
  `(#+lla-pinned-vector-cffi with-vector-output-general%
    #+sbcl with-vector-output-sbcl%                          
    (,output ,pointer ,lla-type ,length)
     ,@body))

(define-with-multiple-bindings with-vector-output)


;;; General pinning code
;;;
;;; We effectively simulate pinning by copying the array to- and from
;;; the memory area.

(defmacro with-copied-elements% ((vector pointer lla-type
                                         &optional (length (gensym* '#:length)))
                                 &body body)
  (check-type pointer symbol)
  (once-only (vector lla-type)
    (with-unique-names (size index)
      `(let ((,length (length ,vector))
             (,size (foreign-size* ,lla-type)))
         (with-foreign-pointer (,pointer (* ,size ,length))
           (dotimes (,index ,length)
             (setf (mem-aref* ,pointer ,lla-type ,index)
                   (coerce* (aref ,vector ,index) ,lla-type)))
           ,@body)))))

(defun copy-elements-from-memory% (pointer lla-type length)
  ;; !! this could be speeded up by conditioning on lla-type
  (let ((vector (lla-vector lla-type length)))
    (dotimes (index length)
      (setf (aref vector index)
            (mem-aref* pointer lla-type index)))
    vector))

(defmacro with-pinned-vector-general% ((vector pointer lla-type
                                               &optional output)
                                       &body body)
  (check-type pointer symbol)
  (once-only (vector lla-type)
    (with-unique-names (length)
      `(with-copied-elements% (,vector ,pointer ,lla-type ,length)
         (multiple-value-prog1
             (progn ,@body)
           ,@(unless (or (not output) (eq output :copy))
               `((setf ,output 
                       (copy-elements-from-memory% ,pointer ,lla-type
                                                   ,length)))))))))

(defmacro with-vector-output-general% ((output pointer lla-type
                                               length) &body body)
  (check-type pointer symbol)
  (once-only (lla-type length)
    `(with-foreign-pointer (,pointer (* (foreign-size* ,lla-type)
                                        ,length))
       (multiple-value-prog1
           (progn ,@body)
         (setf ,output
               (copy-elements-from-memory% ,pointer
                                           ,lla-type ,length))))))

;;; SBCL specific code

#+sbcl
(defmacro pinned-vector-wrapper-sbcl% ((vector pointer) &body body)
  "Pin the vector and bind pointer to its data during body.  This is a
utility function for implementations with pinning, and should not be
used directly in other files, ie it is not part of the interface."
  (check-type pointer symbol)
  (once-only (vector)
    `(sb-sys:with-pinned-objects (,vector)
       (let ((,pointer (sb-sys:vector-sap ,vector)))
	 ,@body))))

#+sbcl
(defmacro with-pinned-vector-sbcl% ((vector pointer lla-type &optional output)
                                    &body body)
  (once-only (lla-type vector)
    (if output
        (let ((vector-copy (gensym* pointer '#:-copy)))
          `(let ((,vector-copy (copy-vector ,vector ,lla-type)))
             (check-type ,lla-type lla-type)
             (multiple-value-prog1 
                 (pinned-vector-wrapper-sbcl% (,vector-copy ,pointer)
                   ,@body)
               ,@(unless (eq output :copy)
                   `((setf ,output ,vector-copy))))))
        `(pinned-vector-wrapper-sbcl% 
             ((if (and (simple-array1? ,vector)
                       (eq (array-lla-type ,vector)
                           ,lla-type))
                  ,vector
                  (copy-vector ,vector ,lla-type))
              ,pointer)
           ,@body))))

#+sbcl
(defmacro with-vector-output-sbcl% ((output pointer lla-type length)
                                    &body body)
  (check-type pointer symbol)
  (once-only (lla-type length)
    (let ((vector (gensym* pointer)))
      `(let ((,vector (lla-vector ,lla-type ,length)))
         (check-type ,lla-type lla-type)
         (multiple-value-prog1
             (pinned-vector-wrapper-sbcl% (,vector ,pointer)
               ,@body)
           (setf ,output ,vector))))))
