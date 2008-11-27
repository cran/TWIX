#include <Rinternals.h>
#include <Rdefines.h>
#include <R.h>
#include <Rmath.h>


double dmax(double *X,int n);
double SUMV_D(double *x, int N);
int SUMV_I(int *x, int N);
void cumsum(double *x, int nx, double *ans);
void rsort_with_x(double *x, double *indx, int n);


