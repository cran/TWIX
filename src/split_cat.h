#include <Rinternals.h>
#include <Rdefines.h>
#include <R.h>
#include <Rmath.h>

SEXP KR_table(SEXP Y, SEXP X);
SEXP my_factor(SEXP sv);
SEXP split_cat( SEXP X, SEXP rsp, SEXP top, SEXP test, SEXP K, SEXP minbuck);
