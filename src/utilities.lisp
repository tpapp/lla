(in-package :lla)

(defun make-symbol* (&rest args)
  "Build a symbol by concatenating each element of ARGS, and intern it
  in LLA.  Elements can be strings or symbols."
  (intern (apply #'concatenate 'string
                 (mapcar (lambda (arg)
                           (etypecase arg
                             (symbol (symbol-name arg))
                             (string arg)))
                         args))
          'lla))

(defun gensym* (&rest args)
  "A version of GENSYM that concatenates args before generating the
symbol.  Also accepts symbols."
  (gensym (apply #'concatenate 'string
                 (mapcar (lambda (arg)
                           (etypecase arg
                             (symbol (symbol-name arg))
                             (string arg)))
                         args))))

(eval-when (:compile-toplevel :load-toplevel :execute)
  (pushnew :muffle-notes cl:*features*))

(defmacro optimize* ((&rest options) &body body)
  `(locally
       (declare (optimize speed)
                #+muffle-notes
                ,@(if (find :muffle-notes options)
                      #+sbcl '((sb-ext:muffle-conditions sb-ext:compiler-note))
                      #-sbcl nil))
     ,@body))
     

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

(defmacro define-with-multiple-bindings (macro)
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

(defun flatten-array (array)
  "Return a 1D displaced version of array."
  (make-array (array-total-size array)
              :element-type (array-element-type array)
              :displaced-to array))
