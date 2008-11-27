#include <Rinternals.h>
#include <Rdefines.h>
#include <R.h>
#include <Rmath.h>

#include "tw_table.h"
#include "Devleaf.h"


static int SUMV(SEXP x)
{
    int i, sum=0;
    for (i = 0; i < LENGTH(x); i++){
		sum += INTEGER(x)[i];
	}
    return sum;
}


static SEXP node_subset( SEXP SubV, SEXP SDATA )
{
    SEXP result,Ldata,Rdata,LDATA,RDATA,Loc_data,class,row_names;
	
	int N = LENGTH(SubV);
	int NLcut = SUMV(SubV);
	int NDATA = LENGTH(SDATA);
    int i,j,r,l,mode;
		
	if(N < 1 || NDATA < 1)
		return(R_NilValue);
		
	PROTECT(LDATA = allocVector(VECSXP,NDATA));
	PROTECT(RDATA = allocVector(VECSXP,NDATA));
		
    for(i=0; i < NDATA; i++){
		Loc_data = VECTOR_ELT(SDATA, i);
		mode = TYPEOF(Loc_data);
		
		if(LENGTH(Loc_data) < 1)
			return(R_NilValue);
		if(NLcut < 1){
			PROTECT(Ldata = allocVector(mode,0));
		}
		else{
			PROTECT(Ldata = allocVector(mode,NLcut));
		}
		if(N-NLcut < 1){
			PROTECT(Rdata = allocVector(mode,0));
		}
		else{
			PROTECT(Rdata = allocVector(mode,N-NLcut));
		}
		switch (TYPEOF(Loc_data)) {
			case LGLSXP:
			case INTSXP:
					for(j=0,l=0,r=0; j < N; j++){
						if(INTEGER(SubV)[j]){
							INTEGER(Ldata)[l]=INTEGER(Loc_data)[j];
							l++;
						}else{
							INTEGER(Rdata)[r]=INTEGER(Loc_data)[j];
							r++;
							}
					}
				break;
			case REALSXP:
				for(j=0,l=0,r=0; j < N; j++){
					if(INTEGER(SubV)[j]){
						REAL(Ldata)[l]=REAL(Loc_data)[j];
						l++;
					}else{
						REAL(Rdata)[r]=REAL(Loc_data)[j];
						r++;
						}
				}
				break;
			case STRSXP:
				for(j=0,l=0,r=0; j < N; j++){
					if(INTEGER(SubV)[j]){
						SET_STRING_ELT(Ldata,l,STRING_ELT(Loc_data,j));
						l++;
					}else{
						SET_STRING_ELT(Rdata,r,STRING_ELT(Loc_data,j));
						r++;
						}
				}
				break;
		}
		if(inherits(Loc_data,"factor")){
			DUPLICATE_ATTRIB(Ldata,Loc_data);
			DUPLICATE_ATTRIB(Rdata,Loc_data);
		}
		SET_VECTOR_ELT(LDATA,i,Ldata);
		SET_VECTOR_ELT(RDATA,i,Rdata);
		UNPROTECT(2);
	}
	
	namesgets(LDATA,getAttrib(SDATA,R_NamesSymbol));
	namesgets(RDATA,getAttrib(SDATA,R_NamesSymbol));
	
	PROTECT(class = mkString("data.frame"));
    setAttrib(LDATA, R_ClassSymbol, class);
	setAttrib(RDATA, R_ClassSymbol, class);
    UNPROTECT(1);
	
	PROTECT(row_names = allocVector(INTSXP, 2));
	INTEGER(row_names)[0] = NA_INTEGER;
	INTEGER(row_names)[1] = LENGTH(Ldata);
	setAttrib(LDATA, R_RowNamesSymbol, row_names);
	UNPROTECT(1);
	
	PROTECT(row_names = allocVector(INTSXP, 2));
	INTEGER(row_names)[0] = NA_INTEGER;
	INTEGER(row_names)[1] = LENGTH(Rdata);
	setAttrib(RDATA, R_RowNamesSymbol, row_names);
	UNPROTECT(1);
		
    PROTECT(result = allocVector(VECSXP,2));
    SET_VECTOR_ELT(result, 0, LDATA);
    SET_VECTOR_ELT(result, 1, RDATA);
    UNPROTECT(3);
    return(result);  
}




SEXP split_rule(SEXP Splitp, SEXP Splitvar, SEXP DEV, SEXP GDEV, SEXP WHICH, SEXP MINBUCK,
                SEXP MINDEV, SEXP DATA, SEXP TDATA, SEXP GESDATA, SEXP GESTDATA, SEXP LEV )
{
    SEXP result_out,Ldata,Ltdata,ans;
    double dev = REAL(DEV)[0];
    int minbuck = INTEGER(MINBUCK)[0];
    double mindev = REAL(MINDEV)[0];
    int NWh = LENGTH(WHICH);
    int ND = LENGTH(DATA);
    int NTD = LENGTH(TDATA);
    int i,j,k,Lobs=0,Robs=0,sw=0,sob=0;
    

	PROTECT(Ldata = allocVector(INTSXP,ND));
    PROTECT(Ltdata = allocVector(INTSXP,NTD));
    for(i=0; i < ND; i++){
        INTEGER(Ldata)[i]=0;
        }
    for(i=0; i < NTD; i++){
        INTEGER(Ltdata)[i]=0;
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
                        INTEGER(Ldata)[i]=1;
                        Lobs++;
                    }
                    else{
                        Robs++;
                    }
                }
                for(i=0; i < NTD; i++){
                    if(REAL(TDATA)[i] < REAL(WHICH)[0]){
                        INTEGER(Ltdata)[i]=1;
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
                            INTEGER(Ldata)[j]=1;
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
                            INTEGER(Ltdata)[j]=1;
                        }
                    }
                }
                sob = (Lobs >= minbuck) && (Robs >= minbuck);
            }
            else{
                sob=0;
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
		
		rsp = VECTOR_ELT(GESDATA, 0);

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
		UNPROTECT(3);
		setAttrib(Sval, R_NamesSymbol, names_Sval);

		PROTECT(result_out = allocVector(VECSXP,6));
		SET_VECTOR_ELT(result_out, 0, ans);
		SET_VECTOR_ELT(result_out, 1, Ldata);
		SET_VECTOR_ELT(result_out, 2, Ltdata);
		SET_VECTOR_ELT(result_out, 3, node_subset(Ldata, GESDATA));
		SET_VECTOR_ELT(result_out, 4, node_subset(Ltdata, GESTDATA));
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

