#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _ncpen_native_cpp_ncpen_fun_(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ncpen_native_cpp_nr_fun_(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ncpen_native_cpp_obj_fun_(SEXP, SEXP, SEXP, SEXP);
extern SEXP _ncpen_native_cpp_obj_grad_fun_(SEXP, SEXP, SEXP, SEXP);
extern SEXP _ncpen_native_cpp_obj_hess_fun_(SEXP, SEXP, SEXP, SEXP);
extern SEXP _ncpen_native_cpp_p_ncpen_fun_(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ncpen_native_cpp_pen_fun_(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ncpen_native_cpp_pen_grad_fun_(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ncpen_native_cpp_qlasso_fun_(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ncpen_native_cpp_set_dev_mode_(SEXP);

static const R_CallMethodDef CallEntries[] = {
     {"_ncpen_native_cpp_ncpen_fun_",    (DL_FUNC) &_ncpen_native_cpp_ncpen_fun_,    22},
     {"_ncpen_native_cpp_nr_fun_",       (DL_FUNC) &_ncpen_native_cpp_nr_fun_,        5},
     {"_ncpen_native_cpp_obj_fun_",      (DL_FUNC) &_ncpen_native_cpp_obj_fun_,       4},
     {"_ncpen_native_cpp_obj_grad_fun_", (DL_FUNC) &_ncpen_native_cpp_obj_grad_fun_,  4},
     {"_ncpen_native_cpp_obj_hess_fun_", (DL_FUNC) &_ncpen_native_cpp_obj_hess_fun_,  4},
     {"_ncpen_native_cpp_p_ncpen_fun_",  (DL_FUNC) &_ncpen_native_cpp_p_ncpen_fun_,  18},
     {"_ncpen_native_cpp_pen_fun_",      (DL_FUNC) &_ncpen_native_cpp_pen_fun_,       5},
     {"_ncpen_native_cpp_pen_grad_fun_", (DL_FUNC) &_ncpen_native_cpp_pen_grad_fun_,  5},
     {"_ncpen_native_cpp_qlasso_fun_",   (DL_FUNC) &_ncpen_native_cpp_qlasso_fun_,   13},
     {"_ncpen_native_cpp_set_dev_mode_", (DL_FUNC) &_ncpen_native_cpp_set_dev_mode_,  1},
     {NULL, NULL, 0}
};

void R_init_ncpen(DllInfo *dll)
{
     R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
     R_useDynamicSymbols(dll, FALSE);
}
