#include <Rinternals.h>
#include <Rdefines.h>
#include <R.h>
#include <Rmath.h>


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

static double clogn(double x)
 { if (x==0) return(0.0); else return(x*log(x)); }


SEXP split_single(  SEXP sv,SEXP rsp,SEXP NN, SEXP svrks )
{
 SEXP result_out, tcls, cls, which2, dev2, globD;


 double maxD=0.0;
 double D = 0.0;
 int i=0,xx;

 double *v_sv=REAL(sv);
 int *v_rsp=INTEGER(rsp);
 int *v_svrks=INTEGER(svrks);
 int len_svrks=LENGTH(svrks);

 xx=cmax(INTEGER(rsp),LENGTH(rsp));

 PROTECT(result_out = allocVector(VECSXP,3));
 PROTECT(tcls = allocVector(INTSXP,xx));
 PROTECT(cls = allocVector(INTSXP,xx));
 double result[LENGTH(sv)], which[LENGTH(sv)];

 for(i=0;i < xx;i++){
 	INTEGER(tcls)[i]=0;
	INTEGER(cls)[i]=0;
	}
 for(i=0;i < len_svrks;i++){
 	result[i]=0.0;
 	which[i]=0.0;
	INTEGER(tcls)[v_rsp[i]-1]+=1;
 }
 D = clogn(LENGTH(sv));

 i =0;
 double lv=0.0, maxDP=0.0, cd =0.0,TD=0.0;
 int isOpt = 0,eq = 0,lct=0;
 while(i < len_svrks){
	lv = v_sv[v_svrks[i]-1];
	if(isOpt) which[eq-1] = (lv+maxDP)/2;
	while(i < len_svrks && lv == v_sv[v_svrks[i]-1]) {
	   if(INTEGER(cls)[v_rsp[v_svrks[i]-1]-1]== 0) {
		INTEGER(cls)[v_rsp[v_svrks[i]-1]-1]=1;
	   } else {
		INTEGER(cls)[v_rsp[v_svrks[i]-1]-1]+=1;
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
		result[eq]=cd;
		maxDP=lv;
		eq++;
	}
 }

 double LL_rsp=LENGTH(rsp);
 for(i=0;i < xx;i++) {
 	TD -= clogn(INTEGER(cls)[i]/LL_rsp)*LL_rsp;
 }
 UNPROTECT(2);
 if((eq-1)>0)
 	i=eq-1;
 else
 	i=0;
 PROTECT(which2 = allocVector(REALSXP,i));
 PROTECT(dev2 = allocVector(REALSXP,i));
 PROTECT(globD = allocVector(REALSXP,1));

 REAL(globD)[0]=TD;
 for(i=0;i<eq-1;i++){
  	REAL(which2)[i]=which[i];
	REAL(dev2)[i]=result[i];
 }
 UNPROTECT(4);
 SET_VECTOR_ELT(result_out, 0, dev2);
 SET_VECTOR_ELT(result_out, 1, which2);
 SET_VECTOR_ELT(result_out, 2, globD);
 return(result_out);
}
