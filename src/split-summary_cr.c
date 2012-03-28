#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>

static double dmax(double *X,int n){
 int i=0;
 double max=0.0;
 for (i=0;i < n; i++) {
 	if (X[i] > max) {
		max=X[i];
		}
	}
 return(max);
}



SEXP split_summary_cr( SEXP BASE, SEXP tol )
{
    int i,j,k,l,sums=0;
    double xx=0.0;
    
    int N = LENGTH(BASE);
    int *ns = Calloc(N, int);
    
    for(i=0; i < N; i++){
        ns[i]=LENGTH(VECTOR_ELT(VECTOR_ELT(BASE,i),0));
        sums+=ns[i];
    }

	double *Ssumm = Calloc(sums, double);
	double *score = Calloc(sums, double);	
	int *Var_id = Calloc(sums, int);
	int *id_Var = Calloc(sums, int);
	int *i_in = Calloc(sums, int);
	    
    k=0;
    for(j=0; j < N; j++){
        SEXP scor = PROTECT(coerceVector(VECTOR_ELT(VECTOR_ELT(BASE,j),3),REALSXP));
        for(i=0; i < ns[j]; i++){
            Ssumm[k]= REAL(VECTOR_ELT(VECTOR_ELT(BASE,j),0))[i];
            score[k]=REAL(scor)[i];
            i_in[k]=k;
            k++;
        }
		UNPROTECT(1);
    }

	SEXP globD = PROTECT(allocVector(REALSXP,1));
    double v=0.0,glob=0.0;
    for(j=0; j < N; j++){
        glob += REAL(VECTOR_ELT(VECTOR_ELT(BASE,j),1))[0];
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
	l=0;
    xx = dmax(Ssumm,k);
	for(i=0; i < k; i++){
		if(Ssumm[i] >= xx*REAL(tol)[0]){
			Var_id[l]=Var_id[i];
			id_Var[l]=id_Var[i];
			score[l]=score[i];
			Ssumm[l]=Ssumm[i];
			score[l]=2.0*(Ssumm[l]/xx)+score[i]/dmax(score,k);
			l++;
			}
		}
	k=l;
	int *i_in2 = Calloc(k, int);
	for(i=0; i < k; i++){
		i_in2[i]=i;
	}
	rsort_with_index(score,i_in2,k);
	SEXP Ssum = PROTECT(allocVector(REALSXP,k));
	SEXP Var_id2 = PROTECT(allocVector(INTSXP,k));
	SEXP id_Var2 = PROTECT(allocVector(INTSXP,k));
	l=0;
	for(i=k-1; i >= 0; i--){
		INTEGER(Var_id2)[l]=Var_id[i_in2[i]];
		INTEGER(id_Var2)[l]=id_Var[i_in2[i]];
		REAL(Ssum)[l]=Ssumm[i_in2[i]];
		l++;
	}
	SEXP which = PROTECT(allocVector(VECSXP,N));
    for(i=0; i < N; i++){
        SET_VECTOR_ELT(which, i, VECTOR_ELT(VECTOR_ELT(BASE,i),2));
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
