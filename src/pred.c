#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>

SEXP getListElementTWIX(SEXP list, const char *str){

	SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
	int i;
     
	for (i = 0; i < length(list); i++)
		if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
			elmt = VECTOR_ELT(list, i);
			break;
		}
	return elmt;
}


SEXP pred_value( SEXP data, SEXP tree, SEXP xvar, SEXP catlev, SEXP meth){

	SEXP rest,rest_spvar,rest_spoint,names,left,right;
    SEXP clev,Dt;
    int i=0,j=0,rule=0,id;
    int lxvar = LENGTH(xvar);
	int meth_pred = INTEGER(meth)[0];

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
//		PROTECT(Dt = coerceVector(data,INTSXP));
        clev = getListElementTWIX(catlev,CHAR(STRING_ELT(xvar, j)));
        for(i=0; i < LENGTH(rest_spoint); i++){
            if(INTEGER(rest_spoint)[i] == 1){
//              if(strcmp(CHAR(STRING_ELT(names, i)), CHAR(STRING_ELT(clev, INTEGER(Dt)[j]-1))) == 0){
				id = REAL(data)[j];
				if(strcmp(CHAR(STRING_ELT(names, i)), CHAR(STRING_ELT(clev, id-1))) == 0){
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
            return(pred_value(data, left, xvar, catlev, meth));
        }
        else {
			if(meth_pred == 0){
				return(getListElementTWIX(left,"Pred.class"));
			}
			else{
				return(getListElementTWIX(left,"Prob"));
			}
		}
    }
    else {
		if(LENGTH(right) == 3){
            return(pred_value(data, right, xvar, catlev, meth));
        }
        else {
			if(meth_pred == 0){
				return(getListElementTWIX(right,"Pred.class"));
			}
			else{
				return(getListElementTWIX(right,"Prob"));
			}
		}
	}
}



SEXP pred_TWIX( SEXP DATA, SEXP TREE, SEXP XVAR, SEXP cat, SEXP N, SEXP pred_meth, SEXP pred_levels) {
    
    SEXP tree_ans,Stree,result_tree,tmp_ans;
    int NDATA = LENGTH(DATA);//INTEGER(N)[0];
    int ltree = LENGTH(TREE);
	int meth_pred = INTEGER(pred_meth)[0];
	int cat_length = INTEGER(pred_levels)[0];
	if(meth_pred == 0){
		PROTECT(tree_ans = allocMatrix(STRSXP,NDATA,ltree));
		int i,j,k=0;
		for(i=0; i < ltree; i++){
			if(LENGTH(VECTOR_ELT(TREE, i)) < 3){
				Stree = VECTOR_ELT(VECTOR_ELT(TREE, i),0);
				for(j=0; j < NDATA; j++){
					SET_STRING_ELT(tree_ans, k,
						STRING_ELT(pred_value(coerceVector(VECTOR_ELT(DATA,j),REALSXP),Stree,XVAR,cat,pred_meth),0));
					k++;
				}
			}
			else{
				Stree = VECTOR_ELT(TREE, i);
				for(j=0; j < NDATA; j++){
					//SET_STRING_ELT(tree_ans, k, STRING_ELT(VECTOR_ELT(TREE, 5),0));
					SET_STRING_ELT(tree_ans, k, STRING_ELT(getListElementTWIX(Stree, "Pred.class"),0));
					k++;
				}
			}
		}
	}
	else{
		PROTECT(tree_ans = allocVector(VECSXP,ltree));
		int start,i,j,k=0;
		for(i=0; i < ltree; i++){
			if(LENGTH(VECTOR_ELT(TREE, i)) < 3){
				Stree = VECTOR_ELT(VECTOR_ELT(TREE, i),0);
				PROTECT(result_tree = allocMatrix(REALSXP,cat_length,NDATA));
				for(j=0, start=0; j < NDATA; j++, start += cat_length){
					tmp_ans = coerceVector(pred_value(coerceVector(VECTOR_ELT(DATA,j),REALSXP),Stree,XVAR,cat,pred_meth),REALSXP);
					for(k=0; k < cat_length; k++){
						REAL(result_tree)[k+start] = REAL(tmp_ans)[k];
					}
				}
				SET_VECTOR_ELT(tree_ans,i,result_tree);
				UNPROTECT(1);
			}
			else{
				Stree = VECTOR_ELT(TREE, i);
				PROTECT(result_tree = allocMatrix(REALSXP,cat_length,NDATA));
				for(j=0, start=0; j < NDATA; j++, start += cat_length){
					tmp_ans = coerceVector(getListElementTWIX(Stree, "Prob"),REALSXP);
					for(k=0; k < cat_length; k++){
						REAL(result_tree)[k+start] = REAL(tmp_ans)[k];
					}
				}
				SET_VECTOR_ELT(tree_ans,i,result_tree);
				UNPROTECT(1);
			}
		}
	}
    UNPROTECT(1);    
    return(tree_ans);
}
