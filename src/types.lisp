(in-package :lla)

;;;; Type synonyms
;;;
;;; These are synonyms for types generally recognized by LLA.  Note that LLA
;;; works fine with any kind of numbers, but may be faster for array element
;;; types it recognized as it does not have to spend time detecting a common
;;; type.  Think of these types as optimization hints when used in MAKE-ARRAY
;;; etc.

(deftype lla-integer ()
  '(signed-byte #-lla::int64 32 #+lla::int64 64))
(deftype lla-single ()
  'single-float)
(deftype lla-double ()
  'double-float)
(deftype lla-complex-single ()
  '(complex single-float))
(deftype lla-complex-double ()
  '(complex double-float))

;;;; Internal types
;;;
;;; LLA uses these constants and the associated functions to represent types
;;; internally.  They are not exported, so we don't worry about name clashes.
;;; The float types are set up so that the common float type can be found
;;; using LOGIOR.

(defconstant +single+ 0)
(defconstant +double+ 1)
(defconstant +complex-single+ 2)
(defconstant +complex-double+ 3)
(defconstant +integer+ 4)

(deftype internal-type ()
  '(integer 0 4))

(define-constant +internal-types+
    (list +single+ +double+ +complex-single+ +complex-double+ +integer+)
  :test #'equal
  :documentation "List of all internal types.")

(define-constant +float-types+
    (list +single+ +double+ +complex-single+ +complex-double+)
  :test #'equal
  :documentation "List of all internal float types.")

(defun lisp-type (internal-type)
  (check-type internal-type internal-type)
  (eswitch (internal-type)
    (+single+ 'lla-single)
    (+double+ 'lla-double)
    (+complex-single+ 'lla-complex-single)
    (+complex-double+ 'lla-complex-double)
    (+integer+ 'lla-integer)))

;;;; Type classification
;;;
;;; In order to decide which of the 4 BLAS/LAPACK functions to call, LLA
;;; attempts to determine the type of objects it is asked to operate on.
;;; Functions below always return an appropriate float type.

(deftype float-type ()
  '(integer 0 3))

(defun number-float-type (number)
  "Return an (internal) float type for a number."
  (declare (optimize speed))
  (typecase number
    (lla-single +single+)
    (lla-double +double+)
    (lla-complex-single +complex-single+)
    (lla-complex-double +complex-double+)
    (complex +complex-double+)
    (real +double+)
    (t (error 'lla-unhandled-type :object number))))

(defun array-float-type (array)
  "Return an (internal) float type for an array.  O(1) when the type can be
detected from the specialized array element type, O(n) otherwise."
  (declare (optimize speed)
           (inline number-float-type))
  (let+ ((element-type (array-element-type array))
         ((&flet is? (type) (type= element-type type))))
    (cond
      ((is? 'lla-single) +single+)
      ((is? 'lla-double) +double+)
      ((is? 'lla-complex-single) +complex-single+)
      ((is? 'lla-complex-double) +complex-double+)
      (t (reduce #'logior (aops:flatten array) :key #'number-float-type)))))

(defun common-float-type (&rest arrays-or-numbers)
  "Return the common (internal) float type for the arguments."
  (reduce #'logior arrays-or-numbers
          :key (lambda (object)
                 (typecase object
                   (number (number-float-type object))
                   (array (array-float-type object))
                   (t (error 'lla-unhandled-type :object object))))))

;;;; Convenience functions for internal types

(defun complex? (internal-type)
  "True iff the internal type is complex."
  (logbitp 1 internal-type))

(defun absolute-square-type (internal-type)
  "Type of (* x (conjugate x))."
  (check-type internal-type float-type)
  (logand 1 internal-type))

(defun epsilon (internal-type)
  "Return the float epsilon for the given internal float type."
  (check-type internal-type float-type)
  (if (logbitp 0 internal-type)
      double-float-epsilon
      single-float-epsilon))
