#include <Rdefines.h>
#include <R.h>
#include <Rmath.h>

#include "tw_table.h"



SEXP Dev_leaf(SEXP y){
	SEXP result = PROTECT(allocVector(REALSXP,1));
	REAL(result)[0] = 0.0;
	
	if(LENGTH(y) < 1){
		UNPROTECT(1);
		return(result);
	}
	SEXP ytab = PROTECT(tw_table(y));
	if(ytab == R_NilValue){
		UNPROTECT(2);
		return(result);
	}
	SEXP x = PROTECT(coerceVector(ytab,REALSXP));
	double TD=0.0;
	int n = LENGTH(y);
	int nm = LENGTH(ytab);

	for(int j=0; j < nm; j++){
		TD = TD + (double) n * clogn(REAL(x)[j] / (double) n);
	}
	if(TD < 0.0)
		REAL(result)[0]= -TD;
	UNPROTECT(3);
	return(result);
}
