#include <Rinternals.h>
#include <Rdefines.h>
#include <R.h>
#include <Rmath.h>


static double clogn(double x)
 { if (x==0) return(0.0); else return(x*log(x)); }

int graycode(int *gray,int x) {
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


SEXP split_cat( SEXP ttot, SEXP TD, SEXP right,
    SEXP numcat, SEXP numclass, SEXP n ,SEXP Wn)
{
 SEXP result_out, Dev, which;

 int i=0,j=0,p=0,k=0,m=0;
 double rtot=REAL(n)[0],ltot=0;
 int Nc = INTEGER(Wn)[0];
 int Nt=INTEGER(numcat)[0];
 int Ncl=INTEGER(numclass)[0];

 int gray[Nt-1],left[Ncl],ind[Nc],tsplit[Nt-1];

 PROTECT(Dev = allocVector(REALSXP,Nc));

 int *D=INTEGER(ttot);
 int *Right=INTEGER(right);
 int ccnt[Nt][Ncl],split[Nc][Nt];

 for (j=0; j < Nc; j++) {
    REAL(Dev)[j]=0.0;
    }
 for(j=0; j < Ncl; j++) {
    left[j]=0.0;
    for( i=0; i<Nt; i++) {
        ccnt[i][j]= D[p];
        p++;
        }
    }

 for(i=0; i < Nt; i++) {
    gray[i]=1;
    tsplit[i]=1;
    }

 while((i=graycode(gray,Nt)) < Nt) {
    if (tsplit[i] == 0) {
        tsplit[i]=1;
        for( j=0; j < Ncl; j++){
            rtot+= ccnt[i][j];
            ltot-= ccnt[i][j];
            Right[j]+= ccnt[i][j];
            left[j]-= ccnt[i][j];
        }
    } else {
        tsplit[i]=0;
        for( j=0; j < Ncl; j++){
            rtot-= ccnt[i][j];
            ltot+= ccnt[i][j];
            Right[j]-= ccnt[i][j];
            left[j]+= ccnt[i][j];
        }
    }
    if (ltot>0  &&  rtot>0) {
    double temp=0.0;
    for(j=0; j < Ncl; j++){
        temp+= -ltot*clogn(left[j]/ltot);
        }
    for(j=0; j < Ncl; j++){
        temp+= -rtot*clogn(Right[j]/rtot);
        }
    REAL(Dev)[k]= REAL(TD)[0]-temp;
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
 rsort_with_index(REAL(Dev),ind,k);
 PROTECT(which = allocVector(INTSXP,Nc*Nt));
 PROTECT(result_out = allocVector(VECSXP,3));
 for(j=0; j < Nc; j++){
    for(p=0; p < Nt; p++){
        INTEGER(which)[m]=split[ind[j]][p];
        m++;
        }
    }
 SET_VECTOR_ELT(result_out, 0, Dev);
 SET_VECTOR_ELT(result_out, 1, TD);
 SET_VECTOR_ELT(result_out, 2, which);
 UNPROTECT(3);
 return(result_out);
 }
 
 
 
 
