#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>

#include "utils.h"


SEXP split_sum( SEXP TBASE, SEXP tol )
{
	SEXP Ssum, Var_id2, id_Var2;
    int i,j,k,l,sums=0;
    double xx=0.0;
    int N = LENGTH(TBASE);
	int *ns = Calloc(N, int);
	
    for(i=0; i < N; i++){
        ns[i] = LENGTH(VECTOR_ELT(VECTOR_ELT(TBASE,i),0));
        sums += ns[i];
    }
	int *Var_id = Calloc(sums, int);
	int *id_Var = Calloc(sums, int);
	int *i_in = Calloc(sums, int);
	double *Ssumm = Calloc(sums, double);

    k=0;
    for(j=0; j < N; j++){
        for(i=0; i < ns[j]; i++){
            Ssumm[k]= REAL(VECTOR_ELT(VECTOR_ELT(TBASE,j),0))[i];
            i_in[k]=k;
            k++;
        }
    }
		
	SEXP globD = PROTECT(allocVector(REALSXP,1));
    REAL(globD)[0] = REAL(VECTOR_ELT(VECTOR_ELT(TBASE,0),1))[0];
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
	l=0;
	xx = dmax(Ssumm,k);
	for(i=0; i < k; i++){
		if(Ssumm[i] >= xx*REAL(tol)[0]){
			Var_id[l]=Var_id[i];
			id_Var[l]=id_Var[i];
			Ssumm[l]=Ssumm[i];
			l++;
			}
		}
		
	k=l;
	int *i_in2 = Calloc(k, int);
	for(i=0; i < k; i++){
		i_in2[i]=i;
	}
	rsort_with_index(Ssumm,i_in2,k);
	Ssum = PROTECT(allocVector(REALSXP,k));
	Var_id2 = PROTECT(allocVector(INTSXP,k));
	id_Var2 = PROTECT(allocVector(INTSXP,k));
	l=0;
	for(i=k-1; i >= 0; i--){
		INTEGER(Var_id2)[l]=Var_id[i_in2[i]];
		INTEGER(id_Var2)[l]=id_Var[i_in2[i]];
		REAL(Ssum)[l]=Ssumm[i];
		l++;
	}
    SEXP which = PROTECT(allocVector(VECSXP,N));
    for(i=0; i < N; i++){
        SET_VECTOR_ELT(which, i, VECTOR_ELT(VECTOR_ELT(TBASE,i),2));
    }
	Free(Var_id); Free(id_Var); Free(Ssumm); Free(ns);
	Free(i_in); Free(i_in2);
    SEXP result_out = PROTECT(allocVector(VECSXP,5));
    SET_VECTOR_ELT(result_out, 0, Ssum);
    SET_VECTOR_ELT(result_out, 1, globD);
    SET_VECTOR_ELT(result_out, 2, which);
    SET_VECTOR_ELT(result_out, 3, id_Var2);
    SET_VECTOR_ELT(result_out, 4, Var_id2);
    UNPROTECT(6);
    return(result_out);  
}
