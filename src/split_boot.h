#include <Rinternals.h>
#include <Rdefines.h>
#include <R.h>
#include <Rmath.h>

int cmax(int *X, int n);
void SampleNoRepl(int k, int n, int *y, int *x);
void dev_split(double *X, int *Y, int nobs, double *dev, double *split, double *TDev);
void tab_num(double *X, int n, double *levels, double *counts, int *k, double *max_count);
SEXP split_boot(SEXP X, SEXP Y, SEXP nboot, SEXP k, SEXP ctopn);
