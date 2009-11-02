#include <Rdefines.h>
#include <R.h>
#include <Rmath.h>

#include "utils.h"


SEXP tw_table(SEXP Data) {

    int NN = LENGTH(Data);
    int i,j,k;

	SEXP C_levels = getAttrib(Data, R_LevelsSymbol);
	int NL = LENGTH(C_levels);
	if(NN < 1 || NL == 0)
		return(R_NilValue);
	
	SEXP C_count = PROTECT(allocVector(INTSXP,NL));
    for (i=1; i <= NL; i++) {
        INTEGER(C_count)[i-1]=0;
        k=0;
        for (j=0; j < NN; j++) {
            if (i == INTEGER(Data)[j]) {
                INTEGER(C_count)[i-1]=++k;
            }
        }
    }
    SET_NAMES(C_count,C_levels);
    UNPROTECT(1);
    return(C_count);
}


SEXP tw_table_prop(SEXP Data) {

    int NN = LENGTH(Data);
    int i,j,k;

	SEXP C_levels = getAttrib(Data, R_LevelsSymbol);
	int NL = LENGTH(C_levels);
	if(NN < 1 || NL == 0)
		return(R_NilValue);
	
	SEXP C_count = PROTECT(allocVector(REALSXP,NL));
    for (i=1; i <= NL; i++) {
        REAL(C_count)[i-1]=0.0;
        k=0.0;
        for (j=0; j < NN; j++) {
            if (i == INTEGER(Data)[j]) {
                REAL(C_count)[i-1]=++k;
            }
        }
		REAL(C_count)[i-1] = REAL(C_count)[i-1] / (double) NN;
    }
    SET_NAMES(C_count,C_levels);
    UNPROTECT(1);
    return(C_count);
}




SEXP my_which(SEXP X, int which)
{
    double s;
    int i;

	SEXP xx = PROTECT(coerceVector(X, REALSXP));
		
    double *r = REAL(xx);
    int n = LENGTH(xx);
    int indx = NA_INTEGER;

    if(which) { /* which.min */
		s = R_PosInf;
		for (i = 0; i < n; i++)
			if ( !ISNAN(r[i]) && (r[i] < s || indx == NA_INTEGER) ) {
			s = r[i]; indx = i;
			}
    } 
	else { /* which.max */
		s = R_NegInf;
		for (i = 0; i < n; i++)
			if ( !ISNAN(r[i]) && (r[i] > s || indx == NA_INTEGER) ) {
			s = r[i]; indx = i;
			}
    }

    i = (indx != NA_INTEGER);
    SEXP ans = PROTECT(allocVector(INTSXP, i ? 1 : 0));
    if (i) {
		INTEGER(ans)[0] = indx + 1;
		if (getAttrib(xx, R_NamesSymbol) != R_NilValue) { /* preserve names */
			SEXP ansname = PROTECT(ScalarString(STRING_ELT(getAttrib(xx, R_NamesSymbol), indx)));
			setAttrib(ans, R_NamesSymbol, ansname);
			UNPROTECT(1);
		}
    }
    UNPROTECT(2);
    return ans;
}



int my_which_raw(SEXP X, int which)
{
    double s, *r;
    int i, n, indx;

	SEXP xx = PROTECT(coerceVector(X, REALSXP));
		
    r = REAL(xx);
    n = LENGTH(xx);
    indx = NA_INTEGER;

    if(which) { /* which.min */
		s = R_PosInf;
		for (i = 0; i < n; i++)
			if ( !ISNAN(r[i]) && (r[i] < s || indx == NA_INTEGER) ) {
			s = r[i]; indx = i;
			}
    } 
	else { /* which.max */
		s = R_NegInf;
		for (i = 0; i < n; i++)
			if ( !ISNAN(r[i]) && (r[i] > s || indx == NA_INTEGER) ) {
			s = r[i]; indx = i;
			}
    }
    UNPROTECT(1);
    return(indx);
}



int potenz(int x, int y)
{
	if(x > 0 && y > 0){
		return(pow(x,y));
	}
	else{
		return(0);
	}
}
