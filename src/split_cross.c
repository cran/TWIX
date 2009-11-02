#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>

static int cmax(int *X,int n){
	int i,max=0;
	for (i=0;i < n; i++) {
		if (X[i] > max) {
			max=X[i];
			}
		}
	return(max);
}

static double clogn(double x)
{
	if (x <= 0.0) 
		return(0.0); 
	else 
		return(x*log(x));
}

// compute k-fold deviance
//
// sv - a numeric vector
// rsp - a integer vector
// NN - a integer vector
// svrks - a integer vetor - sort id's of sv
// s - xgroup id's
// K - a integer - xval 



SEXP split_cross( SEXP sv, SEXP rsp, SEXP NN, SEXP svrks,
			SEXP s, SEXP K, SEXP minbuck)
{
	static double *baseD,*baseW;
	int minbuk = INTEGER(minbuck)[0] - 1;
	int llen_sv=LENGTH(sv);
	int *ss=INTEGER(s);
	int *xobs=INTEGER(NN);
	int xval=INTEGER(K)[0];
	double maxD=0.0,GD=0.0;
	double D = 0.0;
	int k=0,lang=0;
	int i=0,y=0,l=0,xx=0;

	for( y=0; y < xval; y++ ) {
		int i_in[xobs[y]];
		if( xval > 1 ){
			for(l=0; l < LENGTH(s); l++){
				if(ss[l] != y+1)
					i_in[k++] = l;
				}
		} else {
			for(l=0; l < llen_sv; l++){
				i_in[k++]=l; }
		}
		SEXP in_rsp = PROTECT(allocVector(INTSXP,xobs[y]));
		SEXP in_sv = PROTECT(allocVector(REALSXP,xobs[y]));
		SEXP i_sv = PROTECT(allocVector(REALSXP,xobs[y]));

		for(l=0; l < xobs[y] ;l++){
			INTEGER(in_rsp)[l]=INTEGER(rsp)[i_in[l]];
			REAL(i_sv)[l]=REAL(sv)[i_in[l]];
			REAL(in_sv)[l]=REAL(sv)[i_in[l]];
		}
		for(l=0; l < xobs[y] ;l++) i_in[l]=l;
		rsort_with_index(REAL(in_sv),i_in,xobs[y]);
		k=0;
		double *v_sv=REAL(i_sv);
 		int *v_rsp=INTEGER(in_rsp);
 		int *v_svrks=i_in;
 		int len_svrks=LENGTH(in_sv);
		xx=cmax(INTEGER(in_rsp),LENGTH(in_rsp));
	
		int *tcls = Calloc(xx,int);
		int *cls = Calloc(xx,int);
		double *result = Calloc(len_svrks,double);
		double *which = Calloc(len_svrks,double);

 		for(i=0;i < xx;i++){
 			tcls[i]=0;
			cls[i]=0;
			}
		for(i=0;i < len_svrks;i++){
 			result[i]=0.0;
			which[i]=0.0;
			tcls[v_rsp[i]-1]+=1;
			}
 		D = clogn(LENGTH(i_sv));
		i = 0;
		double lv=0.0, maxDP=0.0, cd=0.0, TD=0.0;
		int isOpt = 0, eq = 0, lct=0;

		while(i < len_svrks){
			lv = v_sv[v_svrks[i]];
			if(isOpt) which[eq-1] = (lv+maxDP)/2;
			while(i < len_svrks && lv == v_sv[v_svrks[i]]) {
				if(cls[v_rsp[v_svrks[i]]-1]== 0) {
					cls[v_rsp[v_svrks[i]]-1]=1;
				} else {
					cls[v_rsp[v_svrks[i]]-1]+=1;
					}
				i++;
				lct++;
			}
			cd = D-clogn(lct)-clogn(len_svrks -lct);
			int j=0;
			while(j < xx) {
				cd+=-clogn(tcls[j])+clogn(cls[j])+clogn(tcls[j]-cls[j]);
				j++;
	        }
			isOpt = (cd != maxD) && (i+1 > minbuk) && (i+1 <= len_svrks-minbuk);
			if (isOpt){
				maxD=cd;
				result[eq]=cd;
				maxDP=lv;
				eq++;
			}
		}
		UNPROTECT(3);
		double LL_rsp=LENGTH(in_rsp);
		for(i=0;i < xx; i++) {
			TD -= clogn(cls[i]/LL_rsp)*LL_rsp;
		}
		if(eq > 0){
			baseD = (double*)Realloc( baseD, (lang+eq), double);
			baseW = (double*)Realloc( baseW, (lang+eq), double);
			int j;
			for (j=0; j < eq; j++) {
				baseD[lang+j]=result[j];
				baseW[lang+j]=which[j];
				}
			lang=lang+eq;
			GD+=TD;
		}
		Free(tcls); Free(cls); Free(result); Free(which);
  	}
	if(lang < 1 || GD == 0.0){
		SEXP ans = PROTECT(allocVector(REALSXP,1));
		SEXP ans2 = PROTECT(allocVector(REALSXP,1));
		REAL(ans)[0] = 0.0;
		SEXP result_out = PROTECT(allocVector(VECSXP,4));
		SET_VECTOR_ELT(result_out, 0, ans);
		SET_VECTOR_ELT(result_out, 1, ans);
		REAL(ans2)[0] = 1.0;
		SET_VECTOR_ELT(result_out, 2, ans2);
		SET_VECTOR_ELT(result_out, 3, ans);
		UNPROTECT(3);
		return(result_out);
	}
	int j;
	double dev[lang],wh[lang];
	int seq1[lang],ind[lang];

	for(j=0; j<lang; j++) ind[j]=j;
	rsort_with_index(baseW,ind,lang);
	i=0,j=0,k=0;

	while( i < lang ) {
		wh[k]=baseW[i];
		dev[k]=baseD[ind[i]];
		seq1[k]=1;
		for(j=i+1; j < lang; j++){
			if( wh[k] == baseW[j]){
				dev[k]+= baseD[ind[j]];
				seq1[k] +=1;
			}
     	}
    	i+= seq1[k++];
    }

	SEXP which2 = PROTECT(allocVector(REALSXP,k));
	SEXP dev2 = PROTECT(allocVector(REALSXP,k));
	SEXP seq = PROTECT(allocVector(INTSXP,k));
	SEXP globD = PROTECT(allocVector(REALSXP,1));
	for(i=0;i<k;i++){
		REAL(which2)[i]=wh[i];
		REAL(dev2)[i]=dev[i]/seq1[i];
		INTEGER(seq)[i]=seq1[i];
		}
	REAL(globD)[0]=GD/xval;
	
	SEXP result_out = PROTECT(allocVector(VECSXP,4));
	SET_VECTOR_ELT(result_out, 0, dev2);
	SET_VECTOR_ELT(result_out, 1, which2);
	SET_VECTOR_ELT(result_out, 2, seq);
	SET_VECTOR_ELT(result_out, 3, globD);
	UNPROTECT(5);
	return(result_out);
}

