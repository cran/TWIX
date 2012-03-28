#include <Rinternals.h>
#include <Rdefines.h>
#include <R.h>
#include <Rmath.h>


/* teked from pnmath Package*/
/*							*/

/* Private header file for use during compilation of Mathlib */
#ifndef MATHLIB_PRIVATE_H
#define MATHLIB_PRIVATE_H

#ifndef MATHLIB_STANDALONE
/* Mathlib in R */
# ifdef HAVE_CONFIG_H
#  include <config.h>
# endif
# if defined(HAVE_GLIBC2) && !defined(_BSD_SOURCE)
/* ensure isnan is declared */
#  define _BSD_SOURCE 1
# endif
# if defined(HAVE_GLIBC2) && !defined(_ISOC99_SOURCE)
/* ensure expm1 and log1p are declared */
#  define _ISOC99_SOURCE 1
# endif
#endif

/* need to add LDOUBLE to Rconfig.h and include Rconfig.h in pnmath.h */
#ifdef DODO
#ifdef HAVE_LONG_DOUBLE
#  define LDOUBLE long double
#else
#  define LDOUBLE double
#endif
#else
#  define LDOUBLE long double
#endif
#endif 
/* standalone */



int d2i_round(double Zahl);
int d2i_floor(double Zahl, int Stellen);
double clogn(double x);
int graycode(int *indy,int x);
double dmax(double *X,int n);
double SUMV_D(double *x, int N);
int SUMV_I(int *x, int N);
void cumsum(const double *x, int nx, double *ans);
void rsort_index(double *x, int *indx, int n);
void c_rsort(double *x, int n);
void rsort_with_x(double *x, double *indx, int n);
void rsort_xyz(double *x, double *y, int *indx, int n);
int nrow(SEXP x);
int ncol(SEXP x);
SEXP mysplit(SEXP X);
void loc_maxima( double *x, int N, int *indx );
void dev_seq( int N, int step, int *indx );
SEXP bayes_comb( SEXP pred, SEXP certanty);
