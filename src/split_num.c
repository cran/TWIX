

#include "utils.h"

static int cmax(int *X,int n)
 {
 int i,max=0;
 for (i=0;i < n; i++) {

 	if (X[i] > max) {
		max=X[i];
		}
	}

 return(max);
 }



SEXP split_num( SEXP sv, SEXP rsp, SEXP meth, 
				SEXP step, SEXP topn, SEXP test, SEXP minbuck)
{

	SEXP result_out, which2, dev2, globD;

	double maxD=0.0;
	double D = 0.0;
	int i=0,j,xx;

	int len_svrks = LENGTH(sv);
	int minbuk = INTEGER(minbuck)[0] - 1;
	
	int *v_svrks = Calloc(len_svrks,int);
	double *svrks = Calloc(len_svrks,double);
	
	for(j=0; j < len_svrks; j++){
		v_svrks[j]=j+1;
		svrks[j]=REAL(sv)[j];
	}
	j=0;
	rsort_index(svrks,v_svrks,len_svrks);
	
	Free(svrks);

	double *v_sv=REAL(sv);
	int *v_rsp=INTEGER(rsp);
	xx=cmax(INTEGER(rsp),len_svrks);

	PROTECT(result_out = allocVector(VECSXP,3));
	SEXP names_result;
	PROTECT(names_result = allocVector(STRSXP, 3));
	SET_STRING_ELT(names_result, 0, mkChar("dev"));
	SET_STRING_ELT(names_result, 1, mkChar("globD"));
	SET_STRING_ELT(names_result, 2, mkChar("which"));
	setAttrib(result_out, R_NamesSymbol, names_result);
	UNPROTECT(1);
	
	int *tcls, *cls;
	double *result, *which;
	
	tcls = Calloc(xx,int);
	cls = Calloc(xx,int);
	result = Calloc(len_svrks,double);
	which = Calloc(len_svrks,double);

	for(i=0;i < xx;i++){
		tcls[i]=0;
		cls[i]=0;
	}
	for(i=0;i < len_svrks;i++){
		result[i]=0.0;
		which[i]=0.0;
		tcls[v_rsp[i]-1]+=1;
	}
	D = clogn(len_svrks);

	i = 0;
	double lv=0.0, maxDP=0.0, cd = 0.0,TD=0.0;
	int isOpt = 0,eq = 0,lct=0;
	while(i < len_svrks){
		lv = v_sv[v_svrks[i]-1];
		if(isOpt) which[eq-1] = (lv+maxDP)/2;
		while(i < len_svrks && lv == v_sv[v_svrks[i]-1]){
			if(cls[v_rsp[v_svrks[i]-1]-1]== 0){
				cls[v_rsp[v_svrks[i]-1]-1]=1;
			}
			else{
				cls[v_rsp[v_svrks[i]-1]-1]+=1;
			}
			i++;
			lct++;
		}
		cd = D-clogn(lct)-clogn(len_svrks -lct);
		int j=0;
        while(j < xx){
			cd+=-clogn(tcls[j])+clogn(cls[j])+clogn(tcls[j]-cls[j]);
			j++;
		}
		if(fabs(cd) == FLT_EPSILON || cd < 0.0)
			cd = 0.0;
		isOpt = (cd != maxD) && (i < len_svrks) && (lct >= minbuk) && (len_svrks-lct >= minbuk);
		if(isOpt){
			maxD=cd;
			result[eq]=cd;
			maxDP=lv;
			eq++;
		}
	}

	double LL_rsp=LENGTH(rsp);
	for(i=0; i < xx; i++) {
		TD -= clogn(cls[i]/LL_rsp)*LL_rsp;
	}
/*	if((eq-1) > 0)
		i = eq-1;
	else
		i = 0;
*/
	i=eq;
	Free(v_svrks); Free(tcls); Free(cls);
	//meth=1 -- methode local
	//meth=0 -- methode deviance
	if(i < 1 || TD == 0.0){
		PROTECT(which2 = allocVector(REALSXP,1));
		PROTECT(dev2 = allocVector(REALSXP,1));
		PROTECT(globD = allocVector(REALSXP,1));
		
		REAL(globD)[0] = 0.0;
		REAL(which2)[0] = 0.0;
		REAL(dev2)[0] = 0.0;
		
		Free(result);Free(which);
		SET_VECTOR_ELT(result_out, 0, dev2);
		SET_VECTOR_ELT(result_out, 2, which2);
		SET_VECTOR_ELT(result_out, 1, globD);
		UNPROTECT(4);
		return(result_out);		
	}
	
	
	int NX=i, count=0;
	int st = INTEGER(step)[0];
	int top = INTEGER(topn)[0]+1;
	
	if(LOGICAL(test)[0] && NX > 0){
		PROTECT(which2 = allocVector(REALSXP,NX));
		PROTECT(dev2 = allocVector(REALSXP,NX));
		PROTECT(globD = allocVector(REALSXP,1));
		
		REAL(globD)[0] = TD;
		for(i=0; i < NX; i++){
			REAL(which2)[i] = which[i];
			REAL(dev2)[i] = result[i];
		}
		Free(result);Free(which);
		SET_VECTOR_ELT(result_out, 0, dev2);
		SET_VECTOR_ELT(result_out, 2, which2);
		SET_VECTOR_ELT(result_out, 1, globD);
		UNPROTECT(4);
		return(result_out);
	}
	int *indx;
	indx = Calloc(NX, int);	
	for(j=0; j < NX; j++){
		indx[j]=0;
	}
	if(REAL(meth)[0] == 1.0 && NX > 3){
		loc_maxima(result,NX,indx);
		count++;
	}
	if(REAL(meth)[0] == 0.0 && NX > 3 && st > 1){
		dev_seq(NX,st,indx);
		count++;
	}
	if(count > 0){
		count=0;
		for(j=0; j < NX; j++){
			if(indx[j] > 0){
				result[count] = result[j];
				which[count] = which[j];
				indx[count] = indx[j];
				count++;
			}
		}
	}
	else{
		count=NX;
	}
	Free(indx);
	rsort_with_x(result,which,count);
	if(count < top)
		top = count;

	PROTECT(which2 = allocVector(REALSXP,top));
	PROTECT(dev2 = allocVector(REALSXP,top));
	PROTECT(globD = allocVector(REALSXP,1));
	
	REAL(globD)[0] = TD;
	if(count-1 > 0){
		for(i=count-1, j=0; i > (count-1)-top; i--){
			REAL(which2)[j] = which[i];
			REAL(dev2)[j] = result[i];
			j++;
		}
	}
	else{
		REAL(which2)[0] = which[0];
		REAL(dev2)[0] = result[0];
	}
	Free(result);Free(which);
	SET_VECTOR_ELT(result_out, 0, dev2);
	SET_VECTOR_ELT(result_out, 2, which2);
	SET_VECTOR_ELT(result_out, 1, globD);
	UNPROTECT(4);
	return(result_out);
}



