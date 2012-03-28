

#include "utils.h"

int cmax(int *X, int n)
{
	int i,max=0;
	for (i=0;i < n; i++) {
		if (X[i] > max) {
			max=X[i];
		}
	}
	return(max);
}

/*
Equal probability sampling; without-replacement case
choose k elements from 0 to n-1
*/

void SampleNoRepl(int k, int n, int *y, int *x)
{
    int i, j;
    for (i = 0; i < n; i++)
		x[i] = i;
    for (i = 0; i < k; i++) {
		j = n * unif_rand();
		y[i] = x[j];
		x[j] = x[--n];
    }
}



void dev_split(double *X, int *Y, int nobs, double *dev, double *split, double *TDev)
{

	double maxD=0.0, tmp_dev=0.0;
	double D = 0.0;
	int i=0,j,xx;
	
	int *v_svrks = Calloc(nobs,int);
	double *svrks = Calloc(nobs,double);
	
	for(j=0; j < nobs; j++){
		v_svrks[j] = j+1;
		svrks[j] = X[j];
	}
	rsort_index(svrks,v_svrks,nobs);
	Free(svrks);

	xx = cmax(Y,nobs);		
	int *tcls = Calloc(xx,int);
	int *cls = Calloc(xx,int);

	for(i=0;i < xx;i++){
		tcls[i]=0;
		cls[i]=0;
	}
	for(i=0;i < nobs;i++){
		tcls[Y[i]-1]+=1;
	}
	D = clogn(nobs);

	i = 0;
	double lv=0.0, maxDP=0.0, cd=0.0, TD=0.0;
	int isOpt=0, eq=0, lct=0;
	while(i < nobs){
		lv = X[v_svrks[i]-1];
		if(isOpt) *split = (lv+maxDP)/2;
		while(i < nobs && lv == X[v_svrks[i]-1]){
			if(cls[Y[v_svrks[i]-1]-1]== 0){
				cls[Y[v_svrks[i]-1]-1]=1;
			}
			else{
				cls[Y[v_svrks[i]-1]-1]+=1;
			}
			i++;lct++;
		}
		cd = D-clogn(lct)-clogn(nobs-lct);
		int j=0;
        while(j < xx){
			cd+=-clogn(tcls[j])+clogn(cls[j])+clogn(tcls[j]-cls[j]);
			j++;
		}
		isOpt = cd > maxD;
		if(isOpt){
			maxD=cd;
			tmp_dev = cd;
			maxDP=lv;
			eq++;
		}
	}

	for(i=0; i < xx; i++) {
		TD -= (double) nobs * clogn( (double) cls[i] / (double) nobs);
	}
	*TDev += TD;
	*dev += tmp_dev;
	Free(v_svrks); Free(tcls); Free(cls);
}

/**
    compute a counts of numerical vector
    *\param X a double vector of length n
	*\param n an integer
    *\param levels a double vector
	*\param counts a double vector
*/ 

void tab_num(double *X, int n, double *levels, double *counts, int *k, double *max_count)
{
	int j,m=0;
	c_rsort(X,n);
	levels[m]=X[0];
	for(j=0; j < n; j++){
		if(X[j] > levels[m]){
			m++;
			levels[m]=X[j];
		}
		counts[m]+=1.0;
		if(counts[m] > *max_count)
			*max_count=counts[m];
	}
	*k=m+1;
	rsort_with_x(counts,levels,*k);
}



SEXP split_boot(SEXP X, SEXP Y, SEXP nboot, SEXP k, SEXP ctopn)
{
	int i,j;
	int N = LENGTH(X);
	int NB = REAL(nboot)[0];
	int K = d2i_round((double) N * 0.632);//INTEGER(k)[0];
	int *y_ind = Calloc(K,int);
	int *x_ind = Calloc(N,int);
	double c_dev=0.0, c_gdev=0.0;

	double *xsub = Calloc(K,double);
	int *ysub = Calloc(K,int);
	double *c_split = Calloc(NB,double);

	GetRNGstate();
	for(j=0; j < NB; j++){
		SampleNoRepl(K, N, y_ind, x_ind);
		for(i=0; i < K; i++){
			xsub[i]=REAL(X)[x_ind[i]];
			ysub[i]=INTEGER(Y)[x_ind[i]];
		}
		dev_split(xsub,ysub,K,&c_dev,&c_split[j],&c_gdev);
	}
	PutRNGstate();
	Free(y_ind);Free(x_ind);Free(xsub);Free(ysub);

	c_dev = c_dev / (double) NB;
	double *level = Calloc(NB,double);
	double *count = Calloc(NB,double);
	int new_N = 0;
	double count_max = 0.0;
	int topn_c = INTEGER(ctopn)[0];

	tab_num(c_split, NB, level, count, &new_N, &count_max);
	if(new_N < topn_c)
		topn_c = new_N;
	SEXP GDev = PROTECT(allocVector(REALSXP,1));
	SEXP Split = PROTECT(allocVector(REALSXP,topn_c));
	SEXP Dev = PROTECT(allocVector(REALSXP,topn_c));
	
	for(j=new_N-1,i=0; j >= new_N-topn_c; j--){
		REAL(Dev)[i] = c_dev * (count[j]/count_max);
		REAL(Split)[i] = level[j];
		i++;
	}
	Free(level);Free(count);Free(c_split);
	REAL(GDev)[0] = c_gdev / (double) NB;
	
	SEXP result = PROTECT(allocVector(VECSXP,3));
	SET_VECTOR_ELT(result, 0, Dev);
	SET_VECTOR_ELT(result, 1, GDev);
	SET_VECTOR_ELT(result, 2, Split);
	UNPROTECT(4);
	return(result);
}




