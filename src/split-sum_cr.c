#include <Rinternals.h>
#include <Rdefines.h>
#include <R.h>
#include <Rmath.h>


static double dmax(double *X,int n) {
    int i,max=0.0;
    for (i=0;i < n; i++) {
        if (X[i] > max) {
            max=X[i];
        }
    }
    return(max);
}

SEXP split_sum_cr( SEXP BASE, SEXP NV, SEXP meth, SEXP tol )
{
    SEXP result_out,Ssum,which,Var_id2,id_Var2;
    SEXP globD,dev,scor,data,Ssumm,score;
    int i,j,k,l,sums=0;
    double xx=0;
    
    int N = LENGTH(BASE);
    int ns[N];
    /*PROTECT(ns = allocVector(INTSXP,N));*/
    PROTECT(globD = allocVector(REALSXP,1));
    
    for(i=0;i < N; i++){
        data = VECTOR_ELT(BASE,i);
        ns[i]=LENGTH(VECTOR_ELT(data,0));
        sums+=ns[i];
    }
    
    Ssumm = allocVector(REALSXP,sums);
    PROTECT(score = allocVector(REALSXP,sums));
    int Var_id[sums],id_Var[sums],i_in[sums];
    
    k=0;
    for(j=0; j < N; j++){
        dev = VECTOR_ELT(VECTOR_ELT(BASE,j),0);
        scor = coerceVector(VECTOR_ELT(VECTOR_ELT(BASE,j),3),REALSXP);
        for(i=0; i < ns[j]; i++){
            REAL(Ssumm)[k]= REAL(dev)[i];
            REAL(score)[k]=REAL(scor)[i];
            i_in[k]=k;
            k++;
        }
    }
    double v=0.0,glob=0.0;
    for(j=0; j < N; j++){
        globD =VECTOR_ELT(VECTOR_ELT(BASE,j),1);
        glob+=REAL(globD)[0];
        v++;
    }
    REAL(globD)[0]=glob/v;
    l=0;
    for(i=0; i < N; i++){
        if(ns[i] != 0){
            for(j=1; j <= ns[i]; j++){
                Var_id[l]=j;
                id_Var[l]=i+1;
                l++;
            }
        }
    }
    xx = dmax(REAL(Ssumm),k);
    if(INTEGER(meth)[0] == 0){
        l=0;
        for(i=0; i < k; i++){
            if(REAL(Ssumm)[i] > xx*REAL(tol)[0]){
                Var_id[l]=Var_id[i];
                id_Var[l]=id_Var[i];
                REAL(score)[l]=REAL(score)[i];
                REAL(Ssumm)[l]=REAL(Ssumm)[i];
                REAL(score)[l]=2.0*(REAL(Ssumm)[l]/xx)+REAL(score)[i]/dmax(REAL(score),k);
                l++;
                }
            }
        k=l;
        int i_in2[k];
        for(i=0; i < k; i++){
            i_in2[i]=i;
            }
        rsort_with_index(REAL(score),i_in2,k);
        PROTECT(Ssum = allocVector(REALSXP,k));
        PROTECT(Var_id2 = allocVector(INTSXP,k));
        PROTECT(id_Var2 = allocVector(INTSXP,k));
        l=0;
        for(i=k-1; i >= 0; i--){
            INTEGER(Var_id2)[l]=Var_id[i_in2[i]];
            INTEGER(id_Var2)[l]=id_Var[i_in2[i]];
            REAL(Ssum)[l]=REAL(Ssumm)[i_in2[i]];
            l++;
        }
    }
    else{
        for(i=0; i < k; i++){
            REAL(score)[i]=2.0*(REAL(Ssumm)[i]/xx)+REAL(score)[i]/dmax(REAL(score),k);
        }
        rsort_with_index(REAL(score),i_in,k);
        PROTECT(Ssum = allocVector(REALSXP,k));
        PROTECT(Var_id2 = allocVector(INTSXP,k));
        PROTECT(id_Var2 = allocVector(INTSXP,k));
        l=0;
        for(i=k-1; i >= 0; i--){
            INTEGER(Var_id2)[l]=Var_id[i_in[i]];
            INTEGER(id_Var2)[l]=id_Var[i_in[i]];
            REAL(Ssum)[l]=REAL(Ssumm)[i_in[i]];
            l++;
        }
    }
    PROTECT(which = allocVector(VECSXP,N));
    for(i=0; i < N; i++){
        SET_VECTOR_ELT(which, i, VECTOR_ELT(VECTOR_ELT(BASE,i),2));
    }
    PROTECT(result_out = allocVector(VECSXP,5));
    SET_VECTOR_ELT(result_out, 0, Ssum);
    SET_VECTOR_ELT(result_out, 1, globD);
    SET_VECTOR_ELT(result_out, 2, which);
    SET_VECTOR_ELT(result_out, 3, id_Var2);
    SET_VECTOR_ELT(result_out, 4, Var_id2);
    UNPROTECT(7);
    return(result_out);  
}
