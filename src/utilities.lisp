(in-package :lla)

(defun make-symbol* (&rest args)
  "build a symbol by concatenating each element of ARGS, and intern it
  in the current package.  Elements can be strings or symbols."
  (intern (apply #'concatenate 'string
                 (mapcar (lambda (arg)
                           (etypecase arg
                             (symbol (symbol-name arg))
                             (string arg)))
                         args))))

;;; (make-symbol* "test" "me")        =>   |testme| , :INTERNAL
;;; (make-symbol* "test" 'metoo "me") =>   |testMETOOme| , :INTERNAL
;;; (make-symbol* "TEsT" 'metoo "me") =>   |TEsTMETOOme| , :INTERNAL

(defmacro define-abstract-class (classname super-list &body body)
  "A wrapper for DEFCLASS that lets you define abstract base classes.
   If you try to instantiate an object of this class, a warning is signaled."
  `(progn
     (defclass ,classname ,super-list ,@body)

     ;; Protect against abstract class instantiation.

     ;; We could remove this programmatically later using a
     ;; compile-time constant (or even check the optimization options
     ;; and remove it if SAFETY is set low enough).
     (defmethod initialize-instance :before ((x ,classname) &key)
       (if (eql (type-of x) ',classname)
	   (warn "~A is an abstract base class and not to be instantiated." 
                 (quote ',classname))))))

(defmacro with-multiple-bindings (macro)
  "Define a version of `macro' with multiple arguments, given as a
list.  Application of `macro' will be nested.  The new name is the 
plural of the old one (generated using format)."
  (let ((plural (intern (format nil "~aS" macro))))
    `(defmacro ,plural (bindings &body body)
       ,(format nil "Multiple binding version of ~(~a~)." macro)
       (if bindings
	   `(,',macro ,(car bindings)
		     (,',plural ,(cdr bindings)
			       ,@body))
	   `(progn ,@body)))))

(deftype dimension ()
   "Type for vector/matrix dimensions, basically a nonnegative fixnum."
  '(integer 0 #.most-positive-fixnum))

(defun check-index (index dimension)
  "Error if index is outside dimension."
  ;; ?? should this be a macro, could this signal more information
  (unless (and (<= 0 index) (< index dimension))
    (error "index ~a is outside [0,~a)" index dimension)))

(defun print-length-truncate (dimension)
  "Return values (min dimension *print-length*) and whether the
constraint is binding."
  (if (<= dimension *print-length*)
      (values dimension nil)
      (values *print-length* t)))

(defmacro with-character ((character pointer) &body body)
  "Place character at the memory location of pointer during body."
  (check-type pointer symbol)
  `(with-foreign-object (,pointer :char)
     (setf (mem-ref ,pointer :unsigned-char) (char-code ,character))
     ,@body))

(with-multiple-bindings with-character)
