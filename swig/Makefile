INCLUDE_DIR=/opt/intel/mkl/include
INCLUDE_FILES=$(INCLUDE_DIR)/mkl_cblas.h $(INCLUDE_DIR)/mkl_lapacke.h
LISP_FILENAME=lla-foreign-interface.lisp

all: $(LISP_FILENAME)

$(LISP_FILENAME): lla.i $(INCLUDE_FILES)
	swig2.0 -cffi -I$(INCLUDE_DIR) $<

.PHONY: clean
clean:
	rm $(LISP_FILENAME)
