#include <string.h>

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>


#include "tw_table.h"
#include "Devleaf.h"

 
SEXP KR_table(const SEXP Y, const SEXP X) {
	
	SEXP X_lev = getAttrib(X, R_LevelsSymbol);
	SEXP Y_lev = getAttrib(Y, R_LevelsSymbol);
    	
	int NX = LENGTH(X_lev);
	int NY = LENGTH(Y_lev);
	int N = LENGTH(X);
    int i;	

	SEXP TAB = PROTECT(allocMatrix(INTSXP, NX, NY));
	
	for(i=0; i < NX*NY; i++){
		INTEGER(TAB)[i] = 0;
	}
	for(i=0; i < N; i++){
		INTEGER(TAB)[((INTEGER(Y)[i]-1)*NX+INTEGER(X)[i])-1]+=1;
	}
	SEXP TAB_names = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(TAB_names, 0, X_lev);
    SET_VECTOR_ELT(TAB_names, 1, Y_lev);
    setAttrib(TAB, R_DimNamesSymbol, TAB_names);
    UNPROTECT(2);
    return(TAB);
}



SEXP my_factor(SEXP sv){

	int i=0, j=0, k=0, count=0;
	
	SEXP XX_lev = getAttrib(sv, R_LevelsSymbol);
	int N_XX = LENGTH(XX_lev);
	SEXP id_in = PROTECT(tw_table(sv));
	for(i=0; i < N_XX; i++){
		if(INTEGER(id_in)[i] > 0)
			count++;
	}
	if(count < N_XX){
		SEXP X_lev_new = PROTECT(allocVector(STRSXP, count));
		for(i=0; i < N_XX; i++){
			if(INTEGER(id_in)[i] > 0){
				SET_STRING_ELT(X_lev_new, k++, STRING_ELT(XX_lev,i));
			}
		}
		SEXP C_X = PROTECT(coerceVector(sv,INTSXP));
		k=0;
		for(i=0; i < N_XX; i++){
			if(INTEGER(id_in)[i] < 1){
				for(j=0; j < LENGTH(sv); j++){
					if(INTEGER(C_X)[j] > i-k){
						INTEGER(C_X)[j] = INTEGER(C_X)[j]-1;
					}
				}
			k++;
			}
		}
		setAttrib(C_X, R_LevelsSymbol, X_lev_new);
		setAttrib(C_X, R_ClassSymbol, getAttrib(sv, R_ClassSymbol));
		UNPROTECT(3);
		return(C_X);
	}
	else{
		UNPROTECT(1);
		return(sv);
	}
}



SEXP split_cat( SEXP X, SEXP rsp, SEXP top, SEXP test, SEXP K, SEXP minbuck)
{
	SEXP sv = PROTECT(my_factor(X));
	int i=0, j=0, p=0, k=0, m=0;
	int c_ntop = INTEGER(top)[0];
	double minb = REAL(AS_NUMERIC(minbuck))[0];
	double rtot=LENGTH(rsp),ltot=0.0;
	SEXP GDEV = PROTECT(Dev_leaf(rsp));
	if(REAL(GDEV)[0] > 0.0){	
		int Nt = LENGTH(getAttrib(sv, R_LevelsSymbol));
		int Nc = 2;
		Nc = potenz(2,Nt-1) - 1;
		
		if(Nc >= 1 && Nt > 1){
			SEXP ttot = PROTECT(KR_table(rsp, sv));
			SEXP right = PROTECT(tw_table(rsp));
			int Ncl = LENGTH(right);
			
			//int *gray = Calloc(Nt-1,int);
			int gray[Nt-1];
			//int *left = Calloc(Ncl,int);
			int left[Ncl];
			//int *ind = Calloc(Nc,int);
			int ind[Nc];
			//int *tsplit = Calloc(Nt-1,int);
			int tsplit[Nt-1];
			//double *dev = Calloc(Nc,double);
			double dev[Nc];
			
			int *Right=INTEGER(right);
			int ccnt[Nt][Ncl],split[Nc][Nt];
			
			/*int **ccnt, **split;
			ccnt = Calloc(Nt,int);
			split = Calloc(Nc,int);
			for(j=0; j < Ncl; j++){
				ccnt[j] = Calloc(Ncl,int);
			}
			for(j=0; j < Nt; j++){
				split[j] = Calloc(Nt,int);
			}			
			*/
			for (j=0; j < Nc; j++) {
				dev[j]=0.0;
			}
			for(j=0; j < Ncl; j++) {
				left[j]=0;
				for(i=0; i < Nt; i++) {
					ccnt[i][j]=INTEGER(ttot)[p];
					p++;
				}
			}
			for(i=0; i < Nt; i++) {
				gray[i]=1;
				tsplit[i]=1;
			}
			m=0,k=0;
	
			while((i=graycode(gray,Nt)) < Nt) {
				if (tsplit[i] == 0) {
					tsplit[i]=1;
					for( j=0; j < Ncl; j++){
						rtot+= ccnt[i][j];
						ltot-= ccnt[i][j];
						Right[j]+= ccnt[i][j];
						left[j]-= ccnt[i][j];
					}
				}
				else {
					tsplit[i]=0;
					for( j=0; j < Ncl; j++){
						rtot-= ccnt[i][j];
						ltot+= ccnt[i][j];
						Right[j]-= ccnt[i][j];
						left[j]+= ccnt[i][j];
					}
				}
				if (ltot >= minb && rtot >= minb) {
				//if (ltot >= 1 && rtot >= 1) {
					double temp=0.0;
					for(j=0; j < Ncl; j++){
						temp+= -ltot*clogn(left[j]/ltot);
					}
					for(j=0; j < Ncl; j++){
						temp+= -rtot*clogn(Right[j]/rtot);
					}
					dev[k]= REAL(GDEV)[0]-temp;
				}
				else {
					dev[k]= 0.0;
				}
				ind[k]=k;
				for(j=0; j < Nt; j++){
					split[m][j]=tsplit[j];
				}
				m++;
				k++;
			}
			m=0;
			revsort(dev,ind,k);
			if(INTEGER(test)[0])
				c_ntop = k;
			if(k < c_ntop)
				c_ntop = k;
			SEXP Dev = PROTECT(allocVector(REALSXP,c_ntop));
			SEXP which = PROTECT(allocMatrix(INTSXP, c_ntop, Nt));
			SEXP score = PROTECT(allocVector(INTSXP,c_ntop));
			SEXP result_out = PROTECT(allocVector(VECSXP,4));
			SEXP names_result = PROTECT(allocVector(STRSXP, 4));
			
			SET_STRING_ELT(names_result, 0, mkChar("dev"));
			SET_STRING_ELT(names_result, 1, mkChar("globD"));
			SET_STRING_ELT(names_result, 2, mkChar("which"));
			SET_STRING_ELT(names_result, 3, mkChar("score"));
			setAttrib(result_out, R_NamesSymbol, names_result);
			UNPROTECT(1);

			for (j=0; j < c_ntop; j++) {
				REAL(Dev)[j]=dev[j];
			}
			if(REAL(K)[0] != 0){
				k = d2i_floor(rtot/(rtot*REAL(K)[0]), 2);
			}
			else{
				k = 1;
			}
			for (j=0; j < c_ntop; j++) {
				INTEGER(score)[j]=k;
			}
			
			for(p=0; p < Nt; p++){
				for(j=0; j < c_ntop; j++){
					INTEGER(which)[m]=split[ind[j]][p];
					m++;
					}
				}
			/*Free(dev); Free(tsplit); Free(ind); Free(left); Free(gray);
			for(j=0; j < Ncl; j++){
				Free(ccnt[j]);
			}
			Free(ccnt);
			for(j=0; j < Nt; j++){
				Free(split[j]);
			}
			Free(split);
			*/
			SET_VECTOR_ELT(result_out, 0, Dev);
			SET_VECTOR_ELT(result_out, 1, GDEV);
			SET_VECTOR_ELT(result_out, 2, which);
			SET_VECTOR_ELT(result_out, 3, score);
			UNPROTECT(8);
			return(result_out);
		}
		else{
			SEXP result_out = PROTECT(allocVector(VECSXP,4));
			SEXP names_result = PROTECT(allocVector(STRSXP, 4));
			SET_STRING_ELT(names_result, 0, mkChar("dev"));
			SET_STRING_ELT(names_result, 1, mkChar("globD"));
			SET_STRING_ELT(names_result, 2, mkChar("which"));
			SET_STRING_ELT(names_result, 3, mkChar("score"));
			setAttrib(result_out, R_NamesSymbol, names_result);
			UNPROTECT(1);
			SEXP Dev = PROTECT(allocVector(REALSXP, 1));
			REAL(Dev)[0] = 0.0;
			SET_VECTOR_ELT(result_out, 0, Dev);
			SET_VECTOR_ELT(result_out, 1, Dev);
			SET_VECTOR_ELT(result_out, 2, Dev);
			SET_VECTOR_ELT(result_out, 3, Dev);
			UNPROTECT(4);
			return(result_out);
		}
	}
	else{
		SEXP result_out = PROTECT(allocVector(VECSXP,4));
		SEXP names_result = PROTECT(allocVector(STRSXP, 4));
		SET_STRING_ELT(names_result, 0, mkChar("dev"));
		SET_STRING_ELT(names_result, 1, mkChar("globD"));
		SET_STRING_ELT(names_result, 2, mkChar("which"));
		SET_STRING_ELT(names_result, 3, mkChar("score"));
		setAttrib(result_out, R_NamesSymbol, names_result);
		UNPROTECT(1);
		SEXP Dev = PROTECT(allocVector(REALSXP, 1));
		REAL(Dev)[0] = 0.0;
		SET_VECTOR_ELT(result_out, 0, Dev);
		SET_VECTOR_ELT(result_out, 1, Dev);
		SET_VECTOR_ELT(result_out, 2, Dev);
		SET_VECTOR_ELT(result_out, 3, Dev);
		UNPROTECT(4);
		return(result_out);
	}	
}
 
 
 
 
