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

SEXP split_cross(  SEXP sv,SEXP rsp,SEXP NN, SEXP svrks,
			SEXP s,SEXP K)
{
	SEXP result_out,globD;
	SEXP result,in_rsp,in_sv,i_sv, which, tcls, cls;
	SEXP which2, dev2, seq;

	static double *baseD,*baseW;

	int start=0;
	int llen_sv=LENGTH(sv);
	int *ss=INTEGER(s);
	int *xobs=INTEGER(NN);
	int xval=INTEGER(K)[0];
	double maxD=0.0,GD=0.0;
	double D = 0.0;
	int k=0,lang=0;
	int i=0,y=0,l=0,xx;

	PROTECT(result_out = allocVector(VECSXP,4));

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
		PROTECT(in_rsp = allocVector(INTSXP,xobs[y]));
		PROTECT(in_sv = allocVector(REALSXP,xobs[y]));
		PROTECT(i_sv = allocVector(REALSXP,xobs[y]));

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

		PROTECT(tcls = allocVector(INTSXP,xx));
		PROTECT(cls = allocVector(INTSXP,xx));
		PROTECT(result = allocVector(REALSXP,LENGTH(i_sv)));
		PROTECT(which = allocVector(REALSXP,LENGTH(i_sv)));

 		for(i=0;i < xx;i++){
 			INTEGER(tcls)[i]=0;
			INTEGER(cls)[i]=0;
			}
		for(i=0;i < len_svrks;i++){
 			REAL(result)[i]=0.0;
			REAL(which)[i]=0.0;
			INTEGER(tcls)[v_rsp[i]-1]+=1;
			}
 		D = clogn(LENGTH(i_sv));
		i =0;
		double lv=0.0, maxDP=0.0, cd =0.0,TD=0.0;
		int isOpt = 0,eq = 0,lct=0;

		while(i < len_svrks){
			lv = v_sv[v_svrks[i]];
			if(isOpt) REAL(which)[eq-1] = (lv+maxDP)/2;
			while(i < len_svrks && lv == v_sv[v_svrks[i]]) {
				if(INTEGER(cls)[v_rsp[v_svrks[i]]-1]== 0) {
					INTEGER(cls)[v_rsp[v_svrks[i]]-1]=1;
				} else {
					INTEGER(cls)[v_rsp[v_svrks[i]]-1]+=1;
					}
				i++;
				lct++;
			}
			cd = D-clogn(lct)-clogn(len_svrks -lct);
			int j=0;
			while(j < LENGTH(cls)) {
				cd+=-clogn(INTEGER(tcls)[j])+clogn(INTEGER(cls)[j])+
					clogn(INTEGER(tcls)[j]-INTEGER(cls)[j]);
				j++;
	        }
			isOpt =(cd != maxD);
			if (isOpt){
				maxD=cd;
				REAL(result)[eq]=cd;
				maxDP=lv;
				eq++;
			}
		}
		lang=lang+(eq-1);
		UNPROTECT(5);

		double LL_rsp=LENGTH(in_rsp);
		for(i=0;i < xx;i++) {
			TD -= clogn(INTEGER(cls)[i]/LL_rsp)*LL_rsp;
			}
		baseD= (double*)Realloc( baseD, (lang+3), double);
		baseW= (double*)Realloc( baseW, (lang+3), double);
		int j;
		for (j=0; j < eq-1; j++) {
			baseD[start++]=REAL(result)[j];
			baseW[start-1]=REAL(which)[j];
			}
		GD+=TD;
   		UNPROTECT(2);
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

	PROTECT(which2 = allocVector(REALSXP,k));
	PROTECT(dev2 = allocVector(REALSXP,k));
	PROTECT(seq = allocVector(INTSXP,k));
	PROTECT(globD = allocVector(REALSXP,1));
	for(i=0;i<k;i++){
		REAL(which2)[i]=wh[i];
		REAL(dev2)[i]=dev[i]/seq1[i];
		INTEGER(seq)[i]=seq1[i];
		}
	REAL(globD)[0]=GD/xval;
	SET_VECTOR_ELT(result_out, 0, dev2);
	SET_VECTOR_ELT(result_out, 1, which2);
	SET_VECTOR_ELT(result_out, 2, seq);
	SET_VECTOR_ELT(result_out, 3, globD);
	UNPROTECT(5);
	return(result_out);
}

