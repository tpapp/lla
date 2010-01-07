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

;;;;
;;;;  CL array displacement
;;;;

(defun find-original-array (array)
  "Find the original parent of a displaced array, return this and the
sum of displaced index offsets."
  (let ((sum-of-offsets 0))
    (tagbody
     check-displacement
       (multiple-value-bind (displaced-to displaced-index-offset)
           (array-displacement array)
         (when displaced-to
           (setf array displaced-to)
           (incf sum-of-offsets displaced-index-offset)
           (go check-displacement))))
    (values array sum-of-offsets)))

(defun displace-array (array dimensions index-offset)
  "Make a displaced array from array with the given dimensions and the
index-offset and the same element-type as array.  Tries to displace
from the original array."
  (multiple-value-bind (original-array sum-of-offsets)
      (find-original-array array)
    (make-array dimensions 
                :element-type (array-element-type array)
                :displaced-to original-array
                :displaced-index-offset (+ sum-of-offsets index-offset))))

(defun flatten-array (array)
  "Return a flat (ie rank 1) displaced version of the array."
  (displace-array array (array-total-size array) 0))
