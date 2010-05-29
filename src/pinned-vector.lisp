(in-package :lla)

;;; Pinned vector macros map a numeric vector to a memory area of
;;; value wich specified LLA-TYPE, either for input, input-output or
;;; output, converting the input if required.

;;; NOTE: everything can and will be made to work with other
;;; implementations, but I am concentrating on SBCL for the moment.

#-sbcl (error "You are not using SBCL!  See note in the source ~
  above.")


;;;; All-in-one wrapper macro for numeric vectors.  Just calls one of
;;;; with-nv-input or with-nv-input-output macros below.
;;;; Implementation specific code should just modify the latter
;;;; macros.

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
  `(with-pinned-vector-general% (,vector ,pointer ,lla-type
                                         :output ,output)
    ,@body))

(define-with-multiple-bindings with-pinned-vector)

(defmacro with-vector-output ((output pointer lla-type length)
                              &body body)
  "Allocate a memory area of the given LLA-TYPE and LENGTH, bind its
  address to POINTER, and copy the contents to OUTPUT after BODY."
  `(with-vector-output-general% (,output ,pointer ,lla-type ,length)
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
                                               &key output)
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

;; #+sbcl
;; (defmacro with-pinned-vector ((vector pointer) &body body)
;;   "Pin the vector and bind pointer to its data during body.  This is a
;; utility function for implementations with pinning, and should not be
;; used directly in other files, ie it is not part of the interface."
;;   (check-type pointer symbol)
;;   (once-only (vector)
;;     `(sb-sys:with-pinned-objects (,vector)
;;        (let ((,pointer (sb-sys:vector-sap ,vector)))
;; 	 ,@body))))

;; (defmacro with-pinned-vector-sbcl-input% ((vector pointer lla-type) &body body)
;;   ()
;;   )



;; (defmacro with-nv-input (((numeric-vector-like &optional keyword output) pointer lla-type)
;;                          &body body)
;;   "Make sure that the contents of NUMERIC-VECTOR-LIKE are available at
;; pointer for the duration of BODY.  Optional argument syntax:
;;   (numeric-vector-like) is when body promises not to change the elements,
;;   (numeric-vector-like :copied) copies the elements into a new location,
;;   (numeric-vector-like :output-to output) copies and makes the output
;;   available as a lisp vector of conforming type."
;;   (cond
;;     ((and (symbolp output) (eq keyword :output-to))
;;      ;; (numeric-vector-like :output-to output)
;;      `(with-nv-input-output (,numeric-vector-like ,output ,pointer ,lla-type)
;;         ,@body))
;;     ((and (eq keyword :copied) (null output))
;;      ;; force copy
;;      `(with-nv-input-only (,numeric-vector-like ,pointer ,lla-type t)
;;         ,@body))
;;     ((null keyword)
;;      ;; don't force copy
;;      `(with-nv-input-only (,numeric-vector-like ,pointer ,lla-type nil)
;;         ,@body))
;;     (t (error "Invalid specification."))))

;; (define-with-multiple-bindings with-nv-input)


;; ;;;;
;; ;;;;  Implementation of wrapper macros.
;; ;;;;

;; ;;; The SBCL implementation just uses pinning, and NV-CONVERT for
;; ;;; conversion.


;; #+sbcl
;; (defmacro with-nv-input-only ((numeric-vector-like pointer lla-type copy-p) &body body)
;;   "Makes sure that the contents of numeric-vector-like are available at
;; pointer for the duration of body, with the type LLA-type (converting
;; if necessary).  Unless COPY-P, the body should NOT change the data at
;; the pointer in any way, it is for reading only."
;;   (check-type pointer symbol)
;;   (once-only (numeric-vector-like lla-type)
;;     (with-unique-names (elements)
;;       `(progn
;; 	 (check-type ,numeric-vector-like numeric-vector-like)
;; 	 (check-type ,lla-type lla-type)
;; 	 (let ((,elements (copy-nv-elements% ,numeric-vector-like ,lla-type ,copy-p)))
;; 	   (with-pinned-vector (,elements ,pointer)
;; 	     ,@body))))))

;; #+sbcl
;; (defmacro with-nv-input-output ((numeric-vector-like output pointer lla-type)
;;                                &body body)
;;   "A combination of with-nv-input and -output: input is always copied
;; \(and converted if necessary), is available during body, and the
;; contents end up in a _lisp vector_ of the appropriate type (ie it can
;; be used as elements in MAKE-NV* and MAKE-MATRIX*), assigned to name."
;;   (check-type output symbol)
;;   (check-type pointer symbol)
;;   (once-only (numeric-vector-like lla-type)
;;     `(progn
;;        (check-type ,numeric-vector-like numeric-vector-like)
;;        (check-type ,lla-type lla-type)
;;        (let* ((,output (copy-nv-elements% ,numeric-vector-like ,lla-type t)))
;;         (with-pinned-vector (,output ,pointer)
;;           ,@body)))))


;; #+sbcl
;; (defmacro with-vector-output ((name length pointer lla-type) &body body)
;;   "Allocates a memory area of the given type and length for the
;; duration of body, and makes sure that the contents are assigned to the
;; variable named name at the end as a Lisp array of compatible type.
;; The elements are initialized to zero.  If LENGTH is zero, a length 0
;; vector is created, and pointer is bound to the null pointer."
;;   (check-type name symbol)
;;   (check-type pointer symbol)
;;   (once-only (length lla-type)
;;     `(progn
;;        (check-type ,lla-type lla-type)
;;        (let ((,name (make-nv-elements ,lla-type ,length)))
;;          (if (zerop ,length)
;;              (let ((,pointer (null-pointer)))
;;                ,@body)
;;              (with-pinned-vector (,name ,pointer)
;;                ,@body))))))

;; (define-with-multiple-bindings with-vector-output)
