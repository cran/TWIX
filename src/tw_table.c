#include <Rdefines.h>
#include <R.h>
#include <Rmath.h>


SEXP tw_table(SEXP Data) {

    SEXP count,levels;
    int NN = LENGTH(Data);
    int i,j,k;

	PROTECT(levels = getAttrib(Data, R_LevelsSymbol));
	int NL = LENGTH(levels);
	if(NN < 1 || NL == 0)
		return(R_NilValue);
	
	PROTECT(count = allocVector(INTSXP,NL));

    for (i=1; i <= NL; i++) {
        INTEGER(count)[i-1]=0;
        k=0;
        for (j=0; j < NN; j++) {
            if (i == INTEGER(Data)[j]) {
                INTEGER(count)[i-1]=++k;
            }
        }
    }
    SET_NAMES(count,levels);
    UNPROTECT(2);
    return(count);
}


SEXP tw_table_prop(SEXP Data) {

    SEXP count,levels;
    int NN = LENGTH(Data);
    int i,j,k;

	PROTECT(levels = getAttrib(Data, R_LevelsSymbol));
	int NL = LENGTH(levels);
	if(NN < 1 || NL == 0)
		return(R_NilValue);
	
	PROTECT(count = allocVector(REALSXP,NL));

    for (i=1; i <= NL; i++) {
        REAL(count)[i-1]=0.0;
        k=0.0;
        for (j=0; j < NN; j++) {
            if (i == INTEGER(Data)[j]) {
                REAL(count)[i-1]=++k;
            }
        }
		REAL(count)[i-1] = REAL(count)[i-1]/NN;
    }
    SET_NAMES(count,levels);
    UNPROTECT(2);
    return(count);
}




int nrow(SEXP x) {
    return(INTEGER(getAttrib(x, R_DimSymbol))[0]);
}

int ncol(SEXP x) {
    return(INTEGER(getAttrib(x, R_DimSymbol))[1]);
}

SEXP my_which(SEXP X, int which)
{
    SEXP xx, ans;
    double s, *r;
    int i, n, indx;

	PROTECT(xx = coerceVector(X, REALSXP));
		
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

    i = (indx != NA_INTEGER);
    PROTECT(ans = allocVector(INTSXP, i ? 1 : 0));
    if (i) {
	INTEGER(ans)[0] = indx + 1;
	if (getAttrib(xx, R_NamesSymbol) != R_NilValue) { /* preserve names */
	    SEXP ansnam;
	    PROTECT(ansnam =
		    ScalarString(STRING_ELT(getAttrib(xx, R_NamesSymbol), indx)));
	    setAttrib(ans, R_NamesSymbol, ansnam);
	    UNPROTECT(1);
	}
    }
    UNPROTECT(2);
    return ans;
}

int my_which_raw(SEXP X, int which)
{
    SEXP xx;
    double s, *r;
    int i, n, indx;

	PROTECT(xx = coerceVector(X, REALSXP));
		
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
