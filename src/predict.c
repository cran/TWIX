#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>

SEXP getListElementTWIX(SEXP list, const char *str){
	SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
	int i;
     
	for (i = 0; i < length(list); i++){
		if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
			elmt = VECTOR_ELT(list, i);
			break;
		}
	}
	return elmt;
}


/*
	x is a list
*/

SEXP dev_node(SEXP c_x){
	SEXP names = getAttrib(c_x, R_NamesSymbol);
	int i, i1=0, i2=0;
	for (i = 0; i < length(c_x); i++){
		if(strcmp(CHAR(STRING_ELT(names, i)), "Dev") == 0){
			i1 = i;
		}
		if(strcmp(CHAR(STRING_ELT(names, i)), "dist") == 0){
			i2 = i;
		}
	}
	SET_VECTOR_ELT(c_x, i1, VECTOR_ELT(c_x, i2));
	return c_x;
}



int isElementList(SEXP list, const char *str){

	SEXP names = getAttrib(list, R_NamesSymbol);
	int i,elmt=1;
     
	for (i = 0; i < length(list); i++){
		if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
			elmt = 0;
			break;
		}
	}
	return elmt;
}



SEXP get_tree(SEXP c_tree, SEXP tree_id, SEXP C_I){
	int root = INTEGER(C_I)[0];
	int *id_tree = INTEGER(tree_id);
	SEXP TR = PROTECT(allocVector(VECSXP, 3));
	SEXP names_tree = PROTECT(allocVector(STRSXP, 3));
	SET_STRING_ELT(names_tree, 0, mkChar("split"));
	SET_STRING_ELT(names_tree, 1, mkChar("left"));
	SET_STRING_ELT(names_tree, 2, mkChar("right"));
	setAttrib(TR, R_NamesSymbol, names_tree);
	UNPROTECT(1);
	
	SET_VECTOR_ELT(TR, 0, VECTOR_ELT(VECTOR_ELT(c_tree,0),id_tree[root]-1));
	(INTEGER(C_I)[0])++;
	if(id_tree[INTEGER(C_I)[0]] != 0){
		SET_VECTOR_ELT(TR, 1, get_tree(VECTOR_ELT(VECTOR_ELT(c_tree,1),id_tree[root]-1),tree_id, C_I));
	}
	else{
		if(isElementList(VECTOR_ELT(VECTOR_ELT(c_tree,1),id_tree[root]-1),"split")){
			SET_VECTOR_ELT(TR, 1, VECTOR_ELT(VECTOR_ELT(c_tree,1),id_tree[root]-1));
		}
		else{
			SET_VECTOR_ELT(TR, 1, dev_node(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(c_tree,1),id_tree[root]-1),0),0)));
		}
	}
	(INTEGER(C_I)[0])++;
	if(id_tree[INTEGER(C_I)[0]] != 0){
		SET_VECTOR_ELT(TR, 2, get_tree(VECTOR_ELT(VECTOR_ELT(c_tree,2),id_tree[root]-1),tree_id, C_I));
	}
	else{
		if(isElementList(VECTOR_ELT(VECTOR_ELT(c_tree,2),id_tree[root]-1),"split")){
			SET_VECTOR_ELT(TR, 2, VECTOR_ELT(VECTOR_ELT(c_tree,2),id_tree[root]-1));
		}
		else{
			SET_VECTOR_ELT(TR, 2, dev_node(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(c_tree,2),id_tree[root]-1),0),0)));
		}
	}
	UNPROTECT(1);
	return TR;
}





SEXP pred_value(SEXP data, SEXP tree, SEXP xvar, SEXP catlev, SEXP meth)
{
	int i=0,j=0,rule=0,id;
    int lxvar = LENGTH(xvar);
	int meth_pred = INTEGER(meth)[0];

    SEXP rest = VECTOR_ELT(tree, 0);
	SEXP left = VECTOR_ELT(tree, 1);
    SEXP right = VECTOR_ELT(tree, 2);
    SEXP rest_spvar = VECTOR_ELT(rest, 1);
    SEXP rest_spoint = VECTOR_ELT(rest, 0);
    SEXP names = getAttrib(rest_spoint, R_NamesSymbol);
    
    for(i=0; i<lxvar; i++) {
        if(strcmp(CHAR(STRING_ELT(rest_spvar, 0)), CHAR(STRING_ELT(xvar, i))) == 0){
            j=i;
        }
    }
    if(LENGTH(rest_spoint) > 1){
        SEXP clev = getListElementTWIX(catlev,CHAR(STRING_ELT(xvar, j)));
        for(i=0; i < LENGTH(rest_spoint); i++){
            if(INTEGER(rest_spoint)[i] == 1){
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
    
    SEXP tree_ans, Stree, tmp_ans, tmp_data;
    int NDATA = LENGTH(DATA);
    int ltree = LENGTH(TREE);
	int meth_pred = INTEGER(pred_meth)[0];
	int cat_length = INTEGER(pred_levels)[0];
	if(meth_pred == 0){
		tree_ans = PROTECT(allocMatrix(STRSXP,NDATA,ltree));
		int i,j,k=0;
		for(i=0; i < ltree; i++){
			if(LENGTH(VECTOR_ELT(TREE, i)) < 3){
				Stree = VECTOR_ELT(VECTOR_ELT(TREE, i),0);
				for(j=0; j < NDATA; j++){
					tmp_data = PROTECT(coerceVector(VECTOR_ELT(DATA,j),REALSXP));
					SET_STRING_ELT(tree_ans, k, STRING_ELT(pred_value(tmp_data,Stree,XVAR,cat,pred_meth),0));
					UNPROTECT(1);
					k++;
				}
			}
			else{
				Stree = VECTOR_ELT(TREE, i);
				for(j=0; j < NDATA; j++){
					SET_STRING_ELT(tree_ans, k, STRING_ELT(getListElementTWIX(Stree, "Pred.class"),0));
					k++;
				}
			}
		}
	}
	else{
		tree_ans = PROTECT(allocVector(VECSXP,ltree));
		int start,i,j,k=0;
		for(i=0; i < ltree; i++){
			if(LENGTH(VECTOR_ELT(TREE, i)) < 3){
				Stree = VECTOR_ELT(VECTOR_ELT(TREE, i),0);
				SEXP result_tree = PROTECT(allocMatrix(REALSXP,cat_length,NDATA));
				for(j=0, start=0; j < NDATA; j++, start += cat_length){
					tmp_data = PROTECT(coerceVector(VECTOR_ELT(DATA,j),REALSXP));
					tmp_ans = PROTECT(coerceVector(pred_value(tmp_data,Stree,XVAR,cat,pred_meth),REALSXP));
					for(k=0; k < cat_length; k++){
						REAL(result_tree)[k+start] = REAL(tmp_ans)[k];
					}
					UNPROTECT(2);
				}
				SET_VECTOR_ELT(tree_ans,i,result_tree);
				UNPROTECT(1);
			}
			else{
				Stree = VECTOR_ELT(TREE, i);
				SEXP result_tree = PROTECT(allocMatrix(REALSXP,cat_length,NDATA));
				for(j=0, start=0; j < NDATA; j++, start += cat_length){
					tmp_ans = PROTECT(coerceVector(getListElementTWIX(Stree, "Prob"),REALSXP));
					for(k=0; k < cat_length; k++){
						REAL(result_tree)[k+start] = REAL(tmp_ans)[k];
					}
					UNPROTECT(1);
				}
				SET_VECTOR_ELT(tree_ans,i,result_tree);
				UNPROTECT(1);
			}
		}
	}
    UNPROTECT(1);    
    return(tree_ans);
}
