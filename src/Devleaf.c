#include <Rdefines.h>
#include <R.h>
#include <Rmath.h>

#include "tw_table.h"

static double clogn(double x)
{
	if (x == 0.0) 
		return(0.0); 
	else 
		return(log(x));
}



SEXP Dev_leaf(SEXP y){
	
	SEXP result, x, ytab;
	PROTECT(result = allocVector(REALSXP,1));
	REAL(result)[0]= 0.0;
	
	if(LENGTH(y) < 1){
		UNPROTECT(1);
		return(result);
	}
	
	ytab = tw_table(y);
	if(ytab == R_NilValue){
		UNPROTECT(1);
		return(result);
	}
		
	PROTECT(x = coerceVector(ytab,REALSXP));
	double TD=0.0;
	int n = LENGTH(y);
	int nm = LENGTH(ytab);
	
	for(int j=0; j < nm; j++){
		TD = TD + REAL(x)[j]*clogn(REAL(x)[j]/n);
		}
	if(TD < 0.0)
		REAL(result)[0]= -TD;
	UNPROTECT(2);
	return(result);
}
