(in-package :lla-tests)

(deftestsuite lla-tests () ())

(defun run-lla-tests ()
  "Run all the tests for LLA."
  (run-tests :suite 'lla-tests))

