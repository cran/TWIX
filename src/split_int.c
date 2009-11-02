#include <Rinternals.h>
#include <Rdefines.h>
#include <R.h>
#include <Rmath.h>

#include "utils.h"
#include "split_cat.h"
#include "split_boot.h"
#include "maxstat.h"


SEXP split_int( SEXP sv, SEXP rsp, SEXP meth, 
				SEXP step, SEXP topn, SEXP test, SEXP minbucket)
{

	SEXP C_which, C_dev, C_globD;

	double maxD=0.0;
	double D = 0.0;
	int i=0, j=0;

	int len_svrks = LENGTH(sv);
	int minbuck = INTEGER(minbucket)[0];
	
	int *v_svrks = Calloc(len_svrks,int);
	double *svrks = Calloc(len_svrks,double);
	
	for(j=0; j < len_svrks; j++){
		v_svrks[j]=j+1;
		svrks[j]=REAL(sv)[j];
	}
	rsort_index(svrks,v_svrks,len_svrks);
	Free(svrks);

	double *v_sv=REAL(sv);
	int *v_rsp=INTEGER(rsp);
	int xx = cmax(v_rsp,len_svrks);
	int *tcls = Calloc(xx,int);
	int *cls = Calloc(xx,int);
	double *result = Calloc(len_svrks,double);
	double *which = Calloc(len_svrks,double);
	
	for(i=0; i < xx; i++){
		tcls[i]=0;
		cls[i]=0;
	}
	for(i=0; i < len_svrks; i++){
		result[i]=0.0;
		which[i]=0.0;
		tcls[v_rsp[i]-1]+=1;
	}
	D = clogn(len_svrks);
	i=0;
	double lv=0.0, maxDP=0.0, cd=0.0, TD=0.0;
	int isOpt = 0, eq = 0, lct=0;
	while(i < len_svrks){
		lv = v_sv[v_svrks[i]-1];		
		if(isOpt) which[eq-1] = (lv+maxDP)/2;
		while(i < len_svrks && lv == v_sv[v_svrks[i]-1]){
			if(cls[v_rsp[v_svrks[i]-1]-1] == 0){
				cls[v_rsp[v_svrks[i]-1]-1]=1;
			}
			else{
				cls[v_rsp[v_svrks[i]-1]-1]+=1;
			}
			i++;
			lct++;
		}
		cd = D-clogn(lct)-clogn(len_svrks-lct);		
		int j=0;
        while(j < xx){
			cd+=-clogn(tcls[j])+clogn(cls[j])+clogn(tcls[j]-cls[j]);
			j++;
		}
		if(fabs(cd-0.0) <= FLT_EPSILON || cd < 0.0)
			cd = 0.0;
		isOpt = (cd != maxD) && (i < len_svrks) && (lct >= minbuck) && (len_svrks-lct >= minbuck);
		if(isOpt){
			maxD=cd;
			result[eq]=cd;
			maxDP=lv;
			eq++;
		}
	}

	double LL_rsp = LENGTH(rsp);
	for(i=0; i < xx; i++) {
		TD -= clogn(cls[i]/LL_rsp)*LL_rsp;
	}
	if(fabs(TD-0.0) <= FLT_EPSILON)
		TD = 0.0;
	i = eq;
	Free(v_svrks); Free(tcls); Free(cls);
	
	SEXP result_out = PROTECT(allocVector(VECSXP,3));
	SEXP names_result = PROTECT(allocVector(STRSXP, 3));
	SET_STRING_ELT(names_result, 0, mkChar("dev"));
	SET_STRING_ELT(names_result, 1, mkChar("globD"));
	SET_STRING_ELT(names_result, 2, mkChar("which"));
	setAttrib(result_out, R_NamesSymbol, names_result);
	UNPROTECT(1);
	if(i < 1 || TD == 0.0){
		PROTECT(C_which = allocVector(REALSXP,1));
		PROTECT(C_dev = allocVector(REALSXP,1));
		PROTECT(C_globD = allocVector(REALSXP,1));
		
		REAL(C_globD)[0] = TD;
		REAL(C_which)[0] = 0.0;
		REAL(C_dev)[0] = 0.0;
		
		Free(result);Free(which);
		SET_VECTOR_ELT(result_out, 0, C_dev);
		SET_VECTOR_ELT(result_out, 2, C_which);
		SET_VECTOR_ELT(result_out, 1, C_globD);
		UNPROTECT(4);
		return(result_out);		
	}
	int NX=i, count=0;
	if(LOGICAL(test)[0] && NX > 0){
		PROTECT(C_which = allocVector(REALSXP,NX));
		PROTECT(C_dev = allocVector(REALSXP,NX));
		PROTECT(C_globD = allocVector(REALSXP,1));
		REAL(C_globD)[0] = TD;
		for(i=0; i < NX; i++){
			REAL(C_which)[i] = which[i];
			REAL(C_dev)[i] = result[i];
		}
		Free(result);Free(which);
		SET_VECTOR_ELT(result_out, 0, C_dev);
		SET_VECTOR_ELT(result_out, 2, C_which);
		SET_VECTOR_ELT(result_out, 1, C_globD);
		UNPROTECT(4);
		return(result_out);
	}
	//meth=1 -- methode local
	//meth=0 -- methode deviance
	int st = INTEGER(step)[0];
	int top = INTEGER(topn)[0];
	int *indx = Calloc(NX, int);

	for(j=0; j < NX; j++){
		indx[j]=0;
	}
	if(INTEGER(meth)[0] == 1 && NX > 3){
		loc_maxima(result,NX,indx);
		count++;
	}
	if(INTEGER(meth)[0] == 0 && NX > 3 && st > 1){
		dev_seq(NX,st,indx);
		count++;
	}
	if(count > 0){
		count=0;
		for(j=0; j < NX; j++){
			if(indx[j] > 0){
				result[count] = result[j];
				which[count] = which[j];
				count++;
			}
		}
	}
	else{
		count=NX;
	}
	Free(indx);
	rsort_with_x(result,which,count);
	if(count < top){
		top = count;
	}
	PROTECT(C_which = allocVector(REALSXP,top));
	PROTECT(C_dev = allocVector(REALSXP,top));
	PROTECT(C_globD = allocVector(REALSXP,1));
	REAL(C_globD)[0] = TD;
	if(count-1 > 0){
		for(i=count-1, j=0; i > (count-1)-top; i--){
			REAL(C_which)[j] = which[i];
			REAL(C_dev)[j] = result[i];
			j++;
		}
	}
	else{
		REAL(C_which)[0] = which[0];
		REAL(C_dev)[0] = result[0];
	}
	Free(result);Free(which);
	SET_VECTOR_ELT(result_out, 0, C_dev);
	SET_VECTOR_ELT(result_out, 2, C_which);
	SET_VECTOR_ELT(result_out, 1, C_globD);
	UNPROTECT(4);
	return(result_out);
}


SEXP var_split_dev( const SEXP data, const SEXP meth, const SEXP step, SEXP topn, 
					const SEXP test, const SEXP K, const SEXP minbucket, const SEXP c_lev )
{
	int n, i;
	n = LENGTH(data);
	
	SEXP ans = PROTECT(allocVector(VECSXP, n-1));
	for(i = 1; i < n; i++) {
		if(isFactor(VECTOR_ELT(data, i))){
			SET_VECTOR_ELT(ans, i-1, 
				split_cat(VECTOR_ELT(data, i), VECTOR_ELT(data, 0), topn, test, K, minbucket));
		}
		else{
			SEXP col_data = PROTECT(coerceVector(VECTOR_ELT(data, i),REALSXP));
			if(REAL(K)[0] > 1.0 && INTEGER(c_lev)[0] < 3){
				SET_VECTOR_ELT(ans, i-1, 
					split_boot(col_data, VECTOR_ELT(data, 0), K, K, topn));
			}
			else{
				SET_VECTOR_ELT(ans, i-1, 
					split_int(col_data, VECTOR_ELT(data, 0),
						meth, step, topn, test, minbucket));
			}
			UNPROTECT(1);
		}
	}
	UNPROTECT(1);
	return(ans);
}


SEXP var_split_adj( const SEXP data, const SEXP minprop, const SEXP maxprop, const SEXP test )
{	
	int n, i;
	n = LENGTH(data);
	
	SEXP ans = PROTECT(allocVector(VECSXP, n-1));
	for(i = 1; i < n; i++) {
		SET_VECTOR_ELT(ans, i-1, 
			maxstat(VECTOR_ELT(data, i), VECTOR_ELT(data, 0), 
					minprop, maxprop, test));
	}
	UNPROTECT(1);
	return(ans);
}

