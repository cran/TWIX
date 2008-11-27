#include <Rdefines.h>
#include <R.h>
#include <Rmath.h>


SEXP tw_table(SEXP Data);
SEXP tw_table_prop(SEXP Data);
int nrow(SEXP x);
int ncol(SEXP x);
SEXP my_which(SEXP X, SEXP which);
int my_which_raw(SEXP X, int which);
int potenz(int x, int y);


