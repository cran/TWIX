#include <Rinternals.h>
#include <Rdefines.h>
#include <R.h>
#include <Rmath.h>

SEXP getListElement(SEXP list, char *str)
     {
       SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
       int i;
     
       for (i = 0; i < length(list); i++)
         if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
           elmt = VECTOR_ELT(list, i);
           break;
         }
       return elmt;
}


SEXP pred_value( SEXP data, SEXP tree, SEXP xvar, SEXP catlev)
{
    SEXP rest,rest_spvar,rest_spoint,names,left,right;
    SEXP clev,Dt;
    int i=0,j=0,rule=0;
    int lxvar = LENGTH(xvar);

    rest = VECTOR_ELT(tree, 0);
    left = VECTOR_ELT(tree, 1);
    right = VECTOR_ELT(tree, 2);
    rest_spvar = VECTOR_ELT(VECTOR_ELT(rest, 0), 1);
    rest_spoint = VECTOR_ELT(VECTOR_ELT(rest, 0), 0);
    names = getAttrib(rest_spoint, R_NamesSymbol);
    
    for(i=0; i<lxvar; i++) {
        if(strcmp(CHAR(STRING_ELT(rest_spvar, 0)), CHAR(STRING_ELT(xvar, i))) == 0){
            j=i;
        }
    }
    if(LENGTH(rest_spoint) > 1){
        Dt = coerceVector(data,INTSXP);
        clev = getListElement(catlev,CHAR(STRING_ELT(xvar, j)));
        for(i=0; i < LENGTH(rest_spoint); i++){
            if(INTEGER(rest_spoint)[i] == 1){
                if(strcmp(CHAR(STRING_ELT(names, i)), CHAR(STRING_ELT(clev, INTEGER(Dt)[j]-1))) == 0){
                    rule = 1;
                }
            }
        }
    }
    else {
        rule = REAL(data)[j] < REAL(rest_spoint)[0];
    }
    if(rule){
        if(LENGTH(left) == 3){
            return(pred_value(data, left, xvar, catlev));
        }
        else {
            return(getListElement(left,"Pred.class"));
        }
    }
    else {
        if(LENGTH(right) == 3){
            return(pred_value(data, right, xvar, catlev));
        }
        else {
            return(getListElement(right,"Pred.class"));
        }
    }
}
 
SEXP pred_TWIX( SEXP DATA, SEXP TREE, SEXP XVAR, SEXP cat, SEXP N ) {
    
    SEXP ans,Stree;
    int NDATA = INTEGER(N)[0];
    int ltree = LENGTH(TREE);
    PROTECT(ans = allocMatrix(STRSXP,NDATA,ltree));
    int i,j,k=0;
    
    for(i=0; i < ltree; i++){
        Stree = VECTOR_ELT(TREE, i);
        for(j=0; j < NDATA; j++){
        SET_STRING_ELT(ans, k,STRING_ELT(pred_value(coerceVector(VECTOR_ELT(DATA,j),REALSXP),Stree,XVAR,cat),0));
        k++;
        }
    }
    UNPROTECT(1);    
    return(ans);
}
