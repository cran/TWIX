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

SEXP split_sum( SEXP BASE, SEXP NV, SEXP meth, SEXP tol )
{
    SEXP result_out,Ssum,which,Var_id2,id_Var2;
    SEXP globD,dev,data,Ssumm;
    int i,j,k,l,sums=0;
    double xx=0;
    
    int N = LENGTH(BASE);
    int ns[N];
    PROTECT(globD = allocVector(REALSXP,1));
    
    for(i=0;i < N; i++){
        data = VECTOR_ELT(BASE,i);
        ns[i]=LENGTH(VECTOR_ELT(data,0));
        sums+=ns[i];
    }
    
    Ssumm = allocVector(REALSXP,sums);
    int Var_id[sums],id_Var[sums],i_in[sums];
    
    k=0;
    for(j=0; j < N; j++){
        dev = VECTOR_ELT(VECTOR_ELT(BASE,j),0);
        for(i=0; i < ns[j]; i++){
            REAL(Ssumm)[k]= REAL(dev)[i];
            i_in[k]=k;
            k++;
        }
    }
    REAL(globD)[0]=REAL(VECTOR_ELT(VECTOR_ELT(BASE,0),1))[0];
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
    
    if(INTEGER(meth)[0] == 0){
        l=0;
        xx = dmax(REAL(Ssumm),k);
        for(i=0; i < k; i++){
            if(REAL(Ssumm)[i] > xx*REAL(tol)[0]){
                Var_id[l]=Var_id[i];
                id_Var[l]=id_Var[i];
                REAL(Ssumm)[l]=REAL(Ssumm)[i];
                l++;
                }
            }
        k=l;
        int i_in2[k];
        for(i=0; i < k; i++){
            i_in2[i]=i;
            }
        rsort_with_index(REAL(Ssumm),i_in2,k);
        PROTECT(Ssum = allocVector(REALSXP,k));
        PROTECT(Var_id2 = allocVector(INTSXP,k));
        PROTECT(id_Var2 = allocVector(INTSXP,k));
        l=0;
        for(i=k-1; i >= 0; i--){
            INTEGER(Var_id2)[l]=Var_id[i_in2[i]];
            INTEGER(id_Var2)[l]=id_Var[i_in2[i]];
            REAL(Ssum)[l]=REAL(Ssumm)[i];
            l++;
        }
    }
    else{
        rsort_with_index(REAL(Ssumm),i_in,k);
        PROTECT(Ssum = allocVector(REALSXP,k));
        PROTECT(Var_id2 = allocVector(INTSXP,k));
        PROTECT(id_Var2 = allocVector(INTSXP,k));
        l=0;
        for(i=k-1; i >= 0; i--){
            INTEGER(Var_id2)[l]=Var_id[i_in[i]];
            INTEGER(id_Var2)[l]=id_Var[i_in[i]];
            REAL(Ssum)[l]=REAL(Ssumm)[i];
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
    UNPROTECT(6);
    return(result_out);  
}


SEXP split_rule( SEXP DEV, SEXP GDEV, SEXP WHICH, SEXP MINBUCK,
                SEXP MINDEV, SEXP DATA, SEXP TDATA )
{
    SEXP result_out,Ldata,Ltdata,ans;
    double dev = REAL(DEV)[0];
    int minbuck = INTEGER(MINBUCK)[0];
    double mindev = REAL(MINDEV)[0];
    int NWh = LENGTH(WHICH);
    int ND = LENGTH(DATA);
    int NTD = LENGTH(TDATA);
    int i,j,k,Lobs=0,Robs=0,sw=0,sob=0;
    
    PROTECT(Ldata = allocVector(INTSXP,ND));
    PROTECT(Ltdata = allocVector(INTSXP,NTD));
    for(i=0; i < ND; i++){
        INTEGER(Ldata)[i]=0;
        }
    for(i=0; i < NTD; i++){
        INTEGER(Ltdata)[i]=0;
        }
    if( LENGTH(DEV) != 0) {
        if(NWh == 1) {
            if(dev != 0.0 ){
                sw = (dev != 0.0) && (dev/REAL(GDEV)[0] > mindev);
            }
            else{
                sw=0;
            }
            if(sw){
                for(i=0; i < ND; i++){
                    if(REAL(DATA)[i] < REAL(WHICH)[0]){
                        INTEGER(Ldata)[i]=1;
                        Lobs++;
                    }
                    else{
                        Robs++;
                    }
                }
                for(i=0; i < NTD; i++){
                    if(REAL(TDATA)[i] < REAL(WHICH)[0]){
                        INTEGER(Ltdata)[i]=1;
                    }
                }
                sob = Lobs > minbuck && Robs > minbuck;
            }
            else{
                sob=0;
            }
        }
        else{
            if(dev != 0.0 ){
                sw = (dev != 0.0) && (dev/REAL(GDEV)[0] > mindev);
                }
                else{
                    sw=0;
                }
            if(sw){
                double Lgroup[NWh];
                for(i=0; i < NWh; i++){
                    Lgroup[i]=0.0;
                }
                k=0;
                for(i=0; i < NWh; i++){
                    if(REAL(WHICH)[i] == 1.0){
                        Lgroup[k]=i+1.0;
                        k++;
                    }
                }
                for(i=0; i < k; i++){
                    for(j=0; j < ND; j++){
                        if(REAL(DATA)[j] == Lgroup[i]){
                            INTEGER(Ldata)[j]=1;
                            Lobs++;
                        }
                        else{
                            Robs++;
                        }
                    }   
                }
                for(i=0; i < k; i++){
                    for(j=0; j < NTD; j++){
                        if(REAL(TDATA)[j] == Lgroup[i]){
                            INTEGER(Ltdata)[j]=1;
                        }
                    }
                }
                sob = Lobs > minbuck && Robs > minbuck;
            }
            else{
                sob=0;
            }
        }
    }
    PROTECT(ans = allocVector(INTSXP,1));
    INTEGER(ans)[0]=sob;
    PROTECT(result_out = allocVector(VECSXP,3));
    SET_VECTOR_ELT(result_out, 0, ans);
    SET_VECTOR_ELT(result_out, 1, Ldata);
    SET_VECTOR_ELT(result_out, 2, Ltdata);
    UNPROTECT(4);
    return(result_out);  
}
