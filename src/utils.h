#include <Rinternals.h>
#include <Rdefines.h>
#include <R.h>
#include <Rmath.h>

int d2i_round(double Zahl);
int d2i_floor(double Zahl, int Stellen);
double clogn(double x);
int graycode(int *indy,int x);
double dmax(double *X,int n);
double SUMV_D(double *x, int N);
int SUMV_I(int *x, int N);
void cumsum(double *x, int nx, double *ans);
void rsort_index(double *x, int *indx, int n);
void c_rsort(double *x, int n);
void rsort_with_x(double *x, double *indx, int n);
int nrow(SEXP x);
int ncol(SEXP x);
SEXP mysplit(SEXP X);
void loc_maxima( double *x, int N, int *indx );
void dev_seq( int N, int step, int *indx );
