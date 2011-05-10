#!/usr/local/bin/sbcl --script
(asdf:load-system "lla")
(asdf:load-system "debug-tools")
(sb-ext:save-lisp-and-die #P"/home/tpapp/src/lisp/lla/lla.core")

;;; this script creates a core that has LLA loaded, this is useful when working on
;;; bugs that make the implementation segfault
