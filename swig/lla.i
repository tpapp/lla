%module "lla-foreign-interface"

%insert("lisphead") %{
(in-package :lla)
%}

%typemap(cin) MKL_INT ":int32";
%typemap(cout) MKL_INT ":int32";

%include "mkl_cblas.h";
%include "mkl_lapacke.h";
