#include <Rdefines.h>
#include <R.h>
#include <Rmath.h>

static int cmax(int *X,int n) {
    int i,max=0;
    for (i=0;i < n; i++) {
        if (X[i] > max) {
            max=X[i];
        }
    }
    return(max);
}

SEXP tw_table(SEXP Data, SEXP Lev) {
    SEXP count;
    int NN = LENGTH(Data);
    int NL = LENGTH(Lev);
    int i,j,k,xx=0;
    xx = cmax(INTEGER(Data),NN);    
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
    SET_NAMES(count,Lev);
    UNPROTECT(1);
    return(count);
}
