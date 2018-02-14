#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _ncpen_native_cpp_ncpen_fun_(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ncpen_native_cpp_obj_fun_(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ncpen_native_cpp_obj_grad_fun_(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ncpen_native_cpp_p_ncpen_fun_(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ncpen_native_cpp_qlasso_fun_(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
     {"_ncpen_native_cpp_ncpen_fun_",    (DL_FUNC) &_ncpen_native_cpp_ncpen_fun_,    17},
     {"_ncpen_native_cpp_obj_fun_",      (DL_FUNC) &_ncpen_native_cpp_obj_fun_,       5},
     {"_ncpen_native_cpp_obj_grad_fun_", (DL_FUNC) &_ncpen_native_cpp_obj_grad_fun_,  5},
     {"_ncpen_native_cpp_p_ncpen_fun_",  (DL_FUNC) &_ncpen_native_cpp_p_ncpen_fun_,  14},
     {"_ncpen_native_cpp_qlasso_fun_",   (DL_FUNC) &_ncpen_native_cpp_qlasso_fun_,   10},
     {NULL, NULL, 0}
};

void R_init_ncpen(DllInfo *dll)
{
     R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
     R_useDynamicSymbols(dll, FALSE);
}
