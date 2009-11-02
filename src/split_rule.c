#include <Rinternals.h>
#include <Rdefines.h>
#include <R.h>
#include <Rmath.h>

#include "tw_table.h"
#include "Devleaf.h"
#include "split_cat.h"


int SUMV(SEXP x)
{
    int i, sum=0;
    for (i = 0; i < LENGTH(x); i++){
		sum += INTEGER(x)[i];
	}
    return sum;
}


SEXP my_env(SEXP symbol, SEXP obj, SEXP value, SEXP h, SEXP rho){

	if(LENGTH(obj) > INTEGER(h)[0]){
		SET_VECTOR_ELT(obj, INTEGER(h)[0], value);
	}
	else{
		SEXP yi = PROTECT(allocList(1));
		obj=listAppend(obj,yi);
	}
	UNPROTECT(1);
	return(obj);
}



void append_list(SEXP X, SEXP val, SEXP symbol, SEXP rho)
{
	int i;
	int nx = LENGTH(X);
	
	SEXP ans = PROTECT(allocVector(VECSXP, nx + 1));
	for (i = 0; i < nx; i++){
	    SET_VECTOR_ELT(ans, i, VECTOR_ELT(X, i));
	}
	SET_VECTOR_ELT(ans, nx, val);
	UNPROTECT(1);
	setVar(install(CHAR(symbol)), ans, rho);
}



SEXP node_subset( SEXP SubV, SEXP SDATA )
{
    SEXP result;
	
	int N = LENGTH(SubV);
	int NLcut = SUMV(SubV);
	int NDATA = LENGTH(SDATA);
    int i,j,r,l;
	if(N < 1 || NDATA < 1)
		return(R_NilValue);
	SEXP C_LDATA = PROTECT(allocVector(VECSXP,NDATA));
	SEXP C_RDATA = PROTECT(allocVector(VECSXP,NDATA));
    for(i=0; i < NDATA; i++){
		SEXP C_Loc_data = PROTECT(VECTOR_ELT(SDATA, i));
		if(LENGTH(C_Loc_data) < 1)
			return(R_NilValue);
/*		if(NLcut < 1){
			SEXP C_Ldata = PROTECT(allocVector(mode,0));
		}
		else{
			PROTECT(C_Ldata = allocVector(mode,NLcut));
		}
		if(N-NLcut < 1){
			PROTECT(C_Rdata = allocVector(mode,0));
		}
		else{
			PROTECT(C_Rdata = allocVector(mode,N-NLcut));
		}
*/
		int mode = TYPEOF(C_Loc_data);
		SEXP C_Ldata = PROTECT(allocVector(mode,NLcut));
		SEXP C_Rdata = PROTECT(allocVector(mode,N-NLcut));
		switch (mode) {
			case LGLSXP:
			case INTSXP:
				for(j=0,l=0,r=0; j < N; j++){
					if(INTEGER(SubV)[j]){
						INTEGER(C_Ldata)[l]=INTEGER(C_Loc_data)[j];
						l++;
					}else{
						INTEGER(C_Rdata)[r]=INTEGER(C_Loc_data)[j];
						r++;
					}
				}
				break;
			case REALSXP:
				for(j=0,l=0,r=0; j < N; j++){
					if(INTEGER(SubV)[j]){
						REAL(C_Ldata)[l]=REAL(C_Loc_data)[j];
						l++;
					}else{
						REAL(C_Rdata)[r]=REAL(C_Loc_data)[j];
						r++;
					}
				}
				break;
			case STRSXP:
				for(j=0,l=0,r=0; j < N; j++){
					if(INTEGER(SubV)[j]){
						SET_STRING_ELT(C_Ldata,l,STRING_ELT(C_Loc_data,j));
						l++;
					}else{
						SET_STRING_ELT(C_Rdata,r,STRING_ELT(C_Loc_data,j));
						r++;
					}
				}
				break;
		}
//		if(inherits(C_Loc_data,"factor")){
		if(isFactor(C_Loc_data)){
			DUPLICATE_ATTRIB(C_Ldata,C_Loc_data);
			DUPLICATE_ATTRIB(C_Rdata,C_Loc_data);
		}
		SET_VECTOR_ELT(C_LDATA,i,C_Ldata);
		SET_VECTOR_ELT(C_RDATA,i,C_Rdata);
		UNPROTECT(3);
	}
	
	namesgets(C_LDATA,getAttrib(SDATA,R_NamesSymbol));
	namesgets(C_RDATA,getAttrib(SDATA,R_NamesSymbol));
	
	SEXP C_class = PROTECT(mkString("data.frame"));
    setAttrib(C_LDATA, R_ClassSymbol, C_class);
	setAttrib(C_RDATA, R_ClassSymbol, C_class);
    UNPROTECT(1);
	
	SEXP C_row_names_L = PROTECT(allocVector(INTSXP, 2));
	INTEGER(C_row_names_L)[0] = NA_INTEGER;
	INTEGER(C_row_names_L)[1] = NLcut;//LENGTH(C_Ldata);
	setAttrib(C_LDATA, R_RowNamesSymbol, C_row_names_L);
	UNPROTECT(1);
	
	SEXP C_row_names_R = PROTECT(allocVector(INTSXP, 2));
	INTEGER(C_row_names_R)[0] = NA_INTEGER;
	INTEGER(C_row_names_R)[1] = N-NLcut;//LENGTH(C_Rdata);
	setAttrib(C_RDATA, R_RowNamesSymbol, C_row_names_R);
	UNPROTECT(1);
		
    PROTECT(result = allocVector(VECSXP,2));
    SET_VECTOR_ELT(result, 0, C_LDATA);
    SET_VECTOR_ELT(result, 1, C_RDATA);
    UNPROTECT(3);
    return(result);  
}


SEXP split_rule(SEXP Splitp, SEXP Splitvar, SEXP DEV, SEXP GDEV, SEXP WHICH, SEXP MINBUCK,
                SEXP MINDEV, SEXP DATA, SEXP TDATA, SEXP GESDATA, SEXP GESTDATA, SEXP LEV )
{
    SEXP result_out, C_Ldata, C_Ltdata, ans;
    double dev = REAL(DEV)[0];
    int minbuck = INTEGER(MINBUCK)[0];
    double mindev = REAL(MINDEV)[0];
    int NWh = LENGTH(WHICH);
    int ND = LENGTH(DATA);
    int NTD = LENGTH(TDATA);
    int i,j,k,Lobs=0,Robs=0,sw=0,sob=0;
    

	PROTECT(C_Ldata = allocVector(INTSXP,ND));
    PROTECT(C_Ltdata = allocVector(INTSXP,NTD));
    for(i=0; i < ND; i++){
        INTEGER(C_Ldata)[i]=0;
        }
    for(i=0; i < NTD; i++){
        INTEGER(C_Ltdata)[i]=0;
        }
    if(LENGTH(DEV) > 0) {
        if(NWh == 1) {
            if(dev != 0.0 ){
                sw = (dev != 0.0) && (dev/REAL(GDEV)[0] > mindev);
            }
            else{
                sw=0;
            }
            if(sw){
                for(i=0; i < ND; i++){
                    if(REAL(DATA)[i] < REAL(WHICH)[0]){
                        INTEGER(C_Ldata)[i]=1;
                        Lobs++;
                    }
                    else{
                        Robs++;
                    }
                }
                for(i=0; i < NTD; i++){
                    if(REAL(TDATA)[i] < REAL(WHICH)[0]){
                        INTEGER(C_Ltdata)[i]=1;
                    }
                }
				Robs = ND - Lobs;
                sob = (Lobs >= minbuck) && (Robs >= minbuck);
            }
            else{
                sob=0;
            }
        }
        else{
            if(dev != 0.0 ){
                sw = (dev != 0.0) && (dev/REAL(GDEV)[0] > mindev);
                }
                else{
                    sw=0;
                }
            if(sw){
                double Lgroup[NWh];
                for(i=0; i < NWh; i++){
                    Lgroup[i]=0.0;
                }
                k=0;
                for(i=0; i < NWh; i++){
                    if(REAL(WHICH)[i] == 1.0){
                        Lgroup[k]=i+1.0;
                        k++;
                    }
                }
                for(i=0; i < k; i++){
                    for(j=0; j < ND; j++){
                        if(REAL(DATA)[j] == Lgroup[i]){
                            INTEGER(C_Ldata)[j]=1;
                            Lobs++;
                        }
                        else{
                            Robs++;
                        }
                    }   
                }
				Robs = ND - Lobs;
                for(i=0; i < k; i++){
                    for(j=0; j < NTD; j++){
                        if(REAL(TDATA)[j] == Lgroup[i]){
                            INTEGER(C_Ltdata)[j]=1;
                        }
                    }
                }
                sob = (Lobs >= minbuck) && (Robs >= minbuck);
            }
            else{
                sob = 0;
            }
        }
    }
	PROTECT(ans = allocVector(INTSXP,1));
	INTEGER(ans)[0]=sob;
	
	if(sob){
		SEXP ylev,rsp,ylevt,rspt,Pred_cl,Sval,names_Sval,Obs,Fit_train,Fit_test;
		
		PROTECT(names_Sval = allocVector(STRSXP, 13));
		SET_STRING_ELT(names_Sval, 0, mkChar("Splitp"));
		SET_STRING_ELT(names_Sval, 1, mkChar("Splitvar"));
		SET_STRING_ELT(names_Sval, 2, mkChar("Obs"));
		SET_STRING_ELT(names_Sval, 3, mkChar("Dev"));
		SET_STRING_ELT(names_Sval, 4, mkChar("Dev.test"));
		SET_STRING_ELT(names_Sval, 5, mkChar("fit.tr"));
		SET_STRING_ELT(names_Sval, 6, mkChar("fit.test"));
		SET_STRING_ELT(names_Sval, 7, mkChar("TD"));
		SET_STRING_ELT(names_Sval, 8, mkChar("Pred.class"));
		SET_STRING_ELT(names_Sval, 9, mkChar("Prob"));
		SET_STRING_ELT(names_Sval, 10, mkChar("dist"));
		SET_STRING_ELT(names_Sval, 11, mkChar("ks.t"));
		SET_STRING_ELT(names_Sval, 12, mkChar("dD"));
		
		PROTECT(Sval = allocVector(VECSXP, 13));
		
		SET_VECTOR_ELT(Sval, 0, Splitp);
		SET_VECTOR_ELT(Sval, 1, Splitvar);
		PROTECT(Obs = allocVector(INTSXP,1));
		INTEGER(Obs)[0] = ND;
				
		SET_VECTOR_ELT(Sval, 2, Obs);									/* Obs */		
		SET_VECTOR_ELT(Sval, 3, DEV);									/* Dev */
		
		PROTECT(rsp = VECTOR_ELT(GESDATA, 0));

		ylev = tw_table(rsp);		
		Pred_cl = getAttrib(my_which(ylev,0), R_NamesSymbol);
				
		if(NTD > 1){
			rspt = VECTOR_ELT(GESTDATA, 0);
			ylevt = tw_table(rspt);
			PROTECT(Fit_test = allocVector(INTSXP,1));
			j = 0;
			for(i=0; i < LENGTH(ylevt); i++) { 
				if(strcmp(CHAR(STRING_ELT(Pred_cl,0)), CHAR(STRING_ELT(getAttrib(ylevt, R_NamesSymbol), i))) == 0){
					j=i;
					INTEGER(Fit_test)[0] = INTEGER(ylevt)[i];
				}
			}
			SET_VECTOR_ELT(Sval, 4, Dev_leaf(rspt));					/* Dev.test */
			SET_VECTOR_ELT(Sval, 6, Fit_test);							/* Fit-test */
		}
		else{
			PROTECT(Fit_test = allocVector(INTSXP,1));
			INTEGER(Fit_test)[0] = 0;
			SET_VECTOR_ELT(Sval, 4, Fit_test);							/* Dev.test */
			SET_VECTOR_ELT(Sval, 6, Fit_test);							/* Fit-test */
		}
				
		PROTECT(Fit_train = allocVector(INTSXP,1));
		i = my_which_raw(ylev,0);
		INTEGER(Fit_train)[0] = INTEGER(ylev)[i];
		SET_VECTOR_ELT(Sval, 5, Fit_train);								/* Fit-train */
		SET_VECTOR_ELT(Sval, 7, GDEV);									/* TD */
		SET_VECTOR_ELT(Sval, 8, Pred_cl);								/* Pred.class */
		SET_VECTOR_ELT(Sval, 9, tw_table_prop(rsp));					/* Prob */		
		SET_VECTOR_ELT(Sval, 10, Dev_leaf(rsp));						/* dist */
		SET_VECTOR_ELT(Sval, 11, LEV);									/* ks.t */
		SET_VECTOR_ELT(Sval, 12, DEV);									/* dD   */
		UNPROTECT(4);
		setAttrib(Sval, R_NamesSymbol, names_Sval);

		PROTECT(result_out = allocVector(VECSXP,6));
		SET_VECTOR_ELT(result_out, 0, ans);
		SET_VECTOR_ELT(result_out, 1, C_Ldata);
		SET_VECTOR_ELT(result_out, 2, C_Ltdata);
		SET_VECTOR_ELT(result_out, 3, node_subset(C_Ldata, GESDATA));
		SET_VECTOR_ELT(result_out, 4, node_subset(C_Ltdata, GESTDATA));
		SET_VECTOR_ELT(result_out, 5, Sval);
		UNPROTECT(6);
		return(result_out);
	}
	else{
		PROTECT(result_out = allocVector(VECSXP,1));
		SET_VECTOR_ELT(result_out, 0, ans);
		UNPROTECT(4);
		return(result_out);
	}
}

SEXP split_rule_bag(SEXP spoint, SEXP DEV, SEXP GDEV, SEXP MINBUCK, SEXP MINDEV,
				SEXP ID_VAR, SEXP VAR_ID, SEXP GESDATA, SEXP GESTDATA, SEXP LEV,
				SEXP K, SEXP SV, SEXP LC, SEXP RC, SEXP LCT, SEXP RCT, SEXP rho)
{
	SEXP ans, C_DATA;
	SEXP C_Splitp, C_TDATA;
	int nr, nc, g=0;
	int id_var = INTEGER(ID_VAR)[0];
	int var_id = INTEGER(VAR_ID)[0];

	if(fabs(REAL(DEV)[0]-0.0) <= FLT_EPSILON){
		ans = PROTECT(allocVector(INTSXP,1));
		INTEGER(ans)[0] = 0;
		UNPROTECT(1);
		return(ans);
	}
	if(isMatrix(spoint)){
	    SEXP Xd = getAttrib(spoint, R_DimSymbol);
		nr = INTEGER(Xd)[0];
		nc = INTEGER(Xd)[1];
		PROTECT(C_Splitp = allocVector(INTSXP,nc));
		for(g=0; g < nc; g++){
			INTEGER(C_Splitp)[g] = INTEGER(spoint)[g*nr+var_id];
		}
	}
	else{
		PROTECT(C_Splitp = allocVector(REALSXP,1));
		REAL(C_Splitp)[0] = REAL(spoint)[var_id];
	}
	SEXP C_Splitvar = PROTECT(allocVector(STRSXP,1));
	SET_STRING_ELT(C_Splitvar,0,STRING_ELT(getAttrib(GESDATA, R_NamesSymbol),id_var));
	SEXP DATA_id = PROTECT(VECTOR_ELT(GESDATA,id_var));
	
	if(isMatrix(spoint)){
		PROTECT(C_DATA = my_factor(DATA_id));
		setAttrib(C_Splitp, R_NamesSymbol, getAttrib(C_DATA, R_LevelsSymbol));
		if(GESTDATA != R_NilValue){
			PROTECT(C_TDATA = my_factor(VECTOR_ELT(GESTDATA,id_var)));
		}
		else{
			PROTECT(C_TDATA = allocVector(INTSXP,1));
			INTEGER(C_TDATA)[0] = 0;
		}
	}
	else{
		PROTECT(C_DATA = DATA_id);
		if(GESTDATA != R_NilValue){
			PROTECT(C_TDATA = VECTOR_ELT(GESTDATA,id_var));
		}
		else{
			PROTECT(C_TDATA = allocVector(REALSXP,1));
			REAL(C_TDATA)[0] = 0.0;
		}
	}
//	UNPROTECT(3);

    double dev = REAL(DEV)[0];
    int minbuck = INTEGER(MINBUCK)[0];
    double mindev = REAL(MINDEV)[0];
	int NWh = LENGTH(C_Splitp);
    int ND = LENGTH(C_DATA);
    int NTD = LENGTH(C_TDATA);
    int i,j,k,Lobs=0,Robs=0,sw=0,sob=0;

	SEXP C_Ldata = PROTECT(allocVector(INTSXP,ND));
    SEXP C_Ltdata = PROTECT(allocVector(INTSXP,NTD));
    for(i=0; i < ND; i++){
        INTEGER(C_Ldata)[i]=0;
        }
    for(i=0; i < NTD; i++){
        INTEGER(C_Ltdata)[i]=0;
        }
    if(LENGTH(DEV) > 0) {
        if(NWh == 1) {
            if(fabs(dev - 0.0) > FLT_EPSILON){
                sw = (dev/REAL(GDEV)[0] > mindev);
            }
            else{
                sw=0;
            }
            if(sw){
                for(i=0; i < ND; i++){
                    if(REAL(C_DATA)[i] < REAL(C_Splitp)[0]){
                        INTEGER(C_Ldata)[i]=1;
                        Lobs++;
                    }
                    else{
                        Robs++;
                    }
                }
                for(i=0; i < NTD; i++){
                    if(REAL(C_TDATA)[i] < REAL(C_Splitp)[0]){
                        INTEGER(C_Ltdata)[i]=1;
                    }
                }
				Robs = ND - Lobs;
                sob = (Lobs >= minbuck) && (Robs >= minbuck);
            }
            else{
                sob=0;
            }
        }
        else{
            if(fabs(dev - 0.0) > FLT_EPSILON){
                sw = (dev/REAL(GDEV)[0] > mindev);
                }
                else{
                    sw=0;
                }
            if(sw){
                int Lgroup[NWh];
                for(i=0; i < NWh; i++){
                    Lgroup[i]=0;
                }
                k=0;
                for(i=0; i < NWh; i++){
                    if(INTEGER(C_Splitp)[i] == 1){
                        Lgroup[k]=i+1;
                        k++;
                    }
                }
                for(i=0; i < k; i++){
                    for(j=0; j < ND; j++){
                        if(INTEGER(C_DATA)[j] == Lgroup[i]){
                            INTEGER(C_Ldata)[j]=1;
                            Lobs++;
                        }
                        else{
                            Robs++;
                        }
                    }   
                }
				Robs = ND - Lobs;
                for(i=0; i < k; i++){
                    for(j=0; j < NTD; j++){
                        if(INTEGER(C_TDATA)[j] == Lgroup[i]){
                            INTEGER(C_Ltdata)[j]=1;
                        }
                    }
                }
                sob = (Lobs >= minbuck) && (Robs >= minbuck);
            }
            else{
                sob = 0;
            }
        }
    }
	PROTECT(ans = allocVector(INTSXP,1));
	INTEGER(ans)[0]=sob;
	if(sob){
		SEXP ylevt,rspt;
		
		SEXP names_Sval = PROTECT(allocVector(STRSXP, 13));
		SET_STRING_ELT(names_Sval, 0, mkChar("Splitp"));
		SET_STRING_ELT(names_Sval, 1, mkChar("Splitvar"));
		SET_STRING_ELT(names_Sval, 2, mkChar("Obs"));
		SET_STRING_ELT(names_Sval, 3, mkChar("Dev"));
		SET_STRING_ELT(names_Sval, 4, mkChar("Dev.test"));
		SET_STRING_ELT(names_Sval, 5, mkChar("fit.tr"));
		SET_STRING_ELT(names_Sval, 6, mkChar("fit.test"));
		SET_STRING_ELT(names_Sval, 7, mkChar("TD"));
		SET_STRING_ELT(names_Sval, 8, mkChar("Pred.class"));
		SET_STRING_ELT(names_Sval, 9, mkChar("Prob"));
		SET_STRING_ELT(names_Sval, 10, mkChar("dist"));
		SET_STRING_ELT(names_Sval, 11, mkChar("ks.t"));
		SET_STRING_ELT(names_Sval, 12, mkChar("dD"));
		
		SEXP C_Sval = PROTECT(allocVector(VECSXP, 13));
		SET_VECTOR_ELT(C_Sval, 0, C_Splitp);
		SET_VECTOR_ELT(C_Sval, 1, C_Splitvar);
		SEXP C_Obs = PROTECT(allocVector(INTSXP,1));
		INTEGER(C_Obs)[0] = ND;
				
		SET_VECTOR_ELT(C_Sval, 2, C_Obs);								/* Obs */
		SET_VECTOR_ELT(C_Sval, 3, DEV);									/* Dev */
		
		SEXP C_rsp = PROTECT(VECTOR_ELT(GESDATA, 0));
		SEXP C_ylev = PROTECT(tw_table(C_rsp));		
		SEXP Pred_cl = PROTECT(getAttrib(my_which(C_ylev,0), R_NamesSymbol));
				
		if(GESTDATA != R_NilValue){
			SEXP C_Fit_test = PROTECT(allocVector(INTSXP,1));
			PROTECT(rspt = VECTOR_ELT(GESTDATA, 0));
			ylevt = tw_table(rspt);
			j = 0;
			for(i=0; i < LENGTH(ylevt); i++) {
				if(strcmp(CHAR(STRING_ELT(Pred_cl,0)), CHAR(STRING_ELT(getAttrib(ylevt, R_NamesSymbol), i))) == 0){
					j=i;
					INTEGER(C_Fit_test)[0] = INTEGER(ylevt)[i];
				}
			}			
			SET_VECTOR_ELT(C_Sval, 4, Dev_leaf(rspt));						/* Dev.test */
			SET_VECTOR_ELT(C_Sval, 6, C_Fit_test);							/* Fit-test */
			UNPROTECT(2);
		}
		else{
			SEXP C_Fit_test = PROTECT(allocVector(INTSXP,1));
			INTEGER(C_Fit_test)[0] = 0;
			SET_VECTOR_ELT(C_Sval, 4, C_Fit_test);							/* Dev.test */
			SET_VECTOR_ELT(C_Sval, 6, C_Fit_test);							/* Fit-test */
			UNPROTECT(1);
		}
				
		SEXP C_Fit_train = PROTECT(allocVector(INTSXP,1));
		i = my_which_raw(C_ylev,0);
		INTEGER(C_Fit_train)[0] = INTEGER(C_ylev)[i];
		SET_VECTOR_ELT(C_Sval, 5, C_Fit_train);								/* Fit-train */
		SET_VECTOR_ELT(C_Sval, 7, GDEV);									/* TD */
		SET_VECTOR_ELT(C_Sval, 8, Pred_cl);									/* Pred.class */
		SET_VECTOR_ELT(C_Sval, 9, tw_table_prop(C_rsp));					/* Prob */		
		SET_VECTOR_ELT(C_Sval, 10, Dev_leaf(C_rsp));						/* dist */
		SET_VECTOR_ELT(C_Sval, 11, LEV);									/* ks.t */
		SET_VECTOR_ELT(C_Sval, 12, DEV);									/* dD   */
		UNPROTECT(5);
		setAttrib(C_Sval, R_NamesSymbol, names_Sval);
		
		SEXP C_hn = PROTECT(allocVector(REALSXP,1));
		REAL(C_hn)[0] = REAL(K)[0] + 1.0;
		setVar(install("h_level"), C_hn, rho);
		UNPROTECT(1);
		
		append_list(SV, C_Sval, mkChar("Sval"), rho);
		SEXP Dsub = PROTECT(node_subset(C_Ldata, GESDATA));
		
		append_list(LC, VECTOR_ELT(Dsub,0), mkChar("Lcut"), rho);
		append_list(RC, VECTOR_ELT(Dsub,1), mkChar("Rcut"), rho);
		
		if(GESTDATA != R_NilValue){
			SEXP DsubT = PROTECT(node_subset(C_Ltdata, GESTDATA));
			append_list(LCT, VECTOR_ELT(DsubT,0), mkChar("Ltest"), rho);
			append_list(RCT, VECTOR_ELT(DsubT,1), mkChar("Rtest"), rho);
			UNPROTECT(1);
		}
//		UNPROTECT(8);
		UNPROTECT(11);
		return(ans);
	}
	else{
//		UNPROTECT(5);
		UNPROTECT(8);
		return(ans);
	}
}


SEXP make_node( SEXP topn, const SEXP S, const SEXP MINBUCK, const SEXP MINDEV,
				const SEXP GESDATA, const SEXP GESTDATA, SEXP LEV, SEXP rho)
{
	int i=0;
	int n_top = LENGTH(topn);
	SEXP ans  = PROTECT(allocVector(INTSXP,n_top));
	SEXP C_id_Var_tm = PROTECT(AS_INTEGER(VECTOR_ELT(S,3)));
	SEXP C_Var_id_tm = PROTECT(AS_INTEGER(VECTOR_ELT(S,4)));
	
	for(i=0; i < n_top; i++){
		SEXP C_mn_dev = PROTECT(allocVector(REALSXP,1));
		SEXP C_mn_id_Var = PROTECT(allocVector(INTSXP,1));
		SEXP C_mn_Var_id = PROTECT(allocVector(INTSXP,1));
		REAL(C_mn_dev)[0] = REAL(VECTOR_ELT(S,0))[i];
		INTEGER(C_mn_id_Var)[0] = INTEGER(C_id_Var_tm)[i];
		INTEGER(C_mn_Var_id)[0] = INTEGER(C_Var_id_tm)[i]-1;
		SEXP C_mn_Sval = PROTECT(findVar(install("Sval"),rho));
		SEXP C_mn_Lcut = PROTECT(findVar(install("Lcut"),rho));
		SEXP C_mn_Rcut = PROTECT(findVar(install("Rcut"),rho));
		SEXP C_mn_Ltest = PROTECT(findVar(install("Ltest"),rho));
		SEXP C_mn_Rtest = PROTECT(findVar(install("Rtest"),rho));
		SEXP C_mn_H = PROTECT(findVar(install("h_level"),rho));
		SEXP C_mn_spoint = PROTECT(VECTOR_ELT(VECTOR_ELT(S,2),INTEGER(C_mn_id_Var)[0]-1));

		INTEGER(ans)[i] = INTEGER(split_rule_bag(C_mn_spoint, C_mn_dev, VECTOR_ELT(S,1), MINBUCK, MINDEV,
							C_mn_id_Var, C_mn_Var_id, GESDATA, GESTDATA, LEV, 
							C_mn_H, C_mn_Sval, C_mn_Lcut, C_mn_Rcut, C_mn_Ltest, C_mn_Rtest, rho))[0];

		UNPROTECT(10);
	}
	UNPROTECT(3);
	return(ans);
}

