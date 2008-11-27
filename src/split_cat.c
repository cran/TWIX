#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>

#include "tw_table.h"
#include "Devleaf.h"


static double clogn(double x)
{
	if (x <= 0.0) 
		return(0.0); 
	else 
		return(x*log(x));
}


static int graycode(int *gray,int x) {
    int i;
    int maxc=x;

    for( i=0; i < (maxc-1); i++) {
        if(gray[i] == 1) {
        gray[i]=2;
        return(i);
        }
    else if(gray[i] == 2) gray[i]=1;
    }
    return(maxc);
 }
 
 

SEXP KR_table(SEXP Y, SEXP X) {

    SEXP X_lev, Y_lev, TAB, TAB_names;
	
	PROTECT(X_lev = getAttrib(X, R_LevelsSymbol));
	PROTECT(Y_lev = getAttrib(Y, R_LevelsSymbol));
    	
	int NX = LENGTH(X_lev);
	int NY = LENGTH(Y_lev);
	int N = LENGTH(X);
    int i;	

	PROTECT(TAB = allocMatrix(INTSXP, NX, NY));
	
	for(i=0; i < NX*NY; i++){
		INTEGER(TAB)[i] = 0;
	}
	for(i=0; i < N; i++){
		INTEGER(TAB)[((INTEGER(Y)[i]-1)*NX+INTEGER(X)[i])-1]+=1;
	}
	PROTECT(TAB_names = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(TAB_names, 0, X_lev);
    SET_VECTOR_ELT(TAB_names, 1, Y_lev);
    setAttrib(TAB, R_DimNamesSymbol, TAB_names);
    UNPROTECT(4);
    return(TAB);
}





SEXP split_cat_in( SEXP rsp, SEXP sv)
{
	SEXP result_out, Dev, GDEV, which, ttot, id_in, right, X_lev_new, X_alt;

	int i=0, j=0, p=0, k=0, m=0, count=0;
	double rtot=LENGTH(rsp),ltot=0.0;
	SEXP XX_lev = getAttrib(sv, R_LevelsSymbol);
	int N_XX = LENGTH(XX_lev);

	GDEV = Dev_leaf(rsp);
	id_in = tw_table(sv);


	for (j=0; j < LENGTH(id_in); j++) {
			if(INTEGER(id_in)[j] > 0){
				count+=1;
			}
		}
	if(REAL(GDEV)[0] > 0.0 && count > 1){
		if(count < N_XX){
			PROTECT(X_lev_new = allocVector(STRSXP, count));
			for(i=0; i < N_XX; i++){
				if(INTEGER(id_in)[i] > 0){
					SET_STRING_ELT(X_lev_new, k, STRING_ELT(XX_lev,i));
					k++;
				}
			}	
			UNPROTECT(1);
			PROTECT(X_alt = sv);
			for(i=0; i < N_XX; i++){
				if(INTEGER(id_in)[i] < 1){
					for(j=0; j < LENGTH(sv); j++){
						if(INTEGER(X_alt)[j] > i){
							INTEGER(sv)[j] = INTEGER(sv)[j]-1;
						}
					}
				}
			}
			UNPROTECT(1);
			setAttrib(sv, R_LevelsSymbol, X_lev_new);
		}
		
		int Nt = count;
		int Nc = 2;
		Nc = potenz(2,Nt-1) - 1;
		
		if(Nc >= 1 && Nt > 1){
			PROTECT(ttot = KR_table(rsp, sv));
			right = tw_table(rsp);
			int Ncl = LENGTH(right);
		
			int gray[Nt-1],left[Ncl],ind[Nc],tsplit[Nt-1];
			PROTECT(Dev = allocVector(REALSXP,Nc));

			int *Right=INTEGER(right);
			int ccnt[Nt][Ncl],split[Nc][Nt];

			for (j=0; j < Nc; j++) {
				REAL(Dev)[j]=0.0;
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
	
//		Rprintf("\n::::::::::::::::::::::::   split_cat  9 ::::::::::::::::::\n");
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
				if (ltot > 0.0  &&  rtot > 0.0) {
					double temp=0.0;
					for(j=0; j < Ncl; j++){
						temp+= -ltot*clogn(left[j]/ltot);
					}
					for(j=0; j < Ncl; j++){
						temp+= -rtot*clogn(Right[j]/rtot);
					}
					REAL(Dev)[k]= REAL(GDEV)[0]-temp;
				}
				else {
					REAL(Dev)[k]= 0.0;
				}
				ind[k]=k;
				for(j=0; j < Nt; j++){
					split[m][j]=tsplit[j];
				}
				m++;
				k++;
			}
			m=0;
			revsort(REAL(Dev),ind,k);
		
			PROTECT(which = allocMatrix(INTSXP, Nc, Nt));
			PROTECT(result_out = allocVector(VECSXP,3));

			for(p=0; p < Nt; p++){
				for(j=0; j < Nc; j++){
					INTEGER(which)[m]=split[ind[j]][p];
					m++;
					}
				}
			SET_VECTOR_ELT(result_out, 0, Dev);
			SET_VECTOR_ELT(result_out, 1, GDEV);
			SET_VECTOR_ELT(result_out, 2, which);
			UNPROTECT(4);
			return(result_out);
		}
		else{
			PROTECT(result_out = allocVector(VECSXP,3));
			PROTECT(Dev = allocVector(REALSXP, 1));
			REAL(Dev)[0] = 0.0;
			SET_VECTOR_ELT(result_out, 0, Dev);
			REAL(Dev)[0] = 1.0;
			SET_VECTOR_ELT(result_out, 1, Dev);
			SET_VECTOR_ELT(result_out, 2, Dev);
			UNPROTECT(2);
			return(result_out);
		}
	}
	else{
		PROTECT(result_out = allocVector(VECSXP,3));
		PROTECT(Dev = allocVector(REALSXP, 1));
		REAL(Dev)[0] = 0.0;
		SET_VECTOR_ELT(result_out, 0, Dev);
		REAL(Dev)[0] = 1.0;
		SET_VECTOR_ELT(result_out, 1, Dev);
		SET_VECTOR_ELT(result_out, 2, Dev);
		UNPROTECT(2);
		return(result_out);
	}	
}
 
 
 
 
