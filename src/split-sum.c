#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>

static double dmax(double *X,int n){
 int i=0;
 double max=0.0;
 for (i=0; i < n; i++) {
 	if (X[i] > max) {
		max=X[i];
		}
	}
 return(max);
}



SEXP split_sum( SEXP TBASE, SEXP meth, SEXP tol )
{
    SEXP result_out,Ssum,which,Var_id2,id_Var2;
    SEXP globD,dev,data,Ssumm;
    int i,j,k,l,sums=0;
    double xx=0.0;
    
    int N = LENGTH(TBASE);
    int ns[N];
    PROTECT(globD = allocVector(REALSXP,1));
    
    for(i=0;i < N; i++){
        data = VECTOR_ELT(TBASE,i);
        ns[i]=LENGTH(VECTOR_ELT(data,0));
        sums+=ns[i];
    }
    
    PROTECT(Ssumm = allocVector(REALSXP,sums));
    int Var_id[sums],id_Var[sums],i_in[sums];
    
    k=0;
    for(j=0; j < N; j++){
        dev = VECTOR_ELT(VECTOR_ELT(TBASE,j),0);
        for(i=0; i < ns[j]; i++){
            REAL(Ssumm)[k]= REAL(dev)[i];
            i_in[k]=k;
            k++;
        }
    }
    REAL(globD)[0]=REAL(VECTOR_ELT(VECTOR_ELT(TBASE,0),1))[0];
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
    if(REAL(meth)[0] == 0){
        l=0;
        xx = dmax(REAL(Ssumm),k);
        for(i=0; i < k; i++){
            if(REAL(Ssumm)[i] >= xx*REAL(tol)[0]){
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
/*  Test tol in "single-method" */
        l=0;
        xx = dmax(REAL(Ssumm),k);
        for(i=0; i < k; i++){
            if(REAL(Ssumm)[i] >= xx*REAL(tol)[0]){
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
/*   anstelle alten i_in jetzt i_in2 */
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
    PROTECT(which = allocVector(VECSXP,N));
    for(i=0; i < N; i++){
        SET_VECTOR_ELT(which, i, VECTOR_ELT(VECTOR_ELT(TBASE,i),2));
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
