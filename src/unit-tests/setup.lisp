(in-package :lla-unit-tests)

;; TEST SUITES

(deftestsuite lla-unit-tests () ())


;; EXTERNAL

(defun run-lla-tests ()
  "Run all the tests for LLA."
  (run-tests :suite 'lla-unit-tests))

;;;; run all tests
;;;; (run-lla-tests)
