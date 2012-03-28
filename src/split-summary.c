
#include "utils.h"



static SEXP cummin(SEXP x, SEXP s)
{
    int i;
    double min, *rx = REAL(x), *rs = REAL(s);
    min = R_PosInf; /* always positive, not NA */
    for (i = 0 ; i < length(x) ; i++ ) {
		if (ISNAN(rx[i]) || ISNAN(min))
			min = min + rx[i];  /* propagate NA and NaN */
		else
			min = (min < rx[i]) ? min : rx[i];
		rs[i] = min;
    }
    return s;
}



static void cummin_C(double *rx, double *rs, int N)
{
    int i;
    double min;
    min = R_PosInf; /* always positive, not NA */
    for (i = 0 ; i < N ; i++ ) {
		if (ISNAN(rx[i]) || ISNAN(min))
			min = min + rx[i];  /* propagate NA and NaN */
		else
			min = (min < rx[i]) ? min : rx[i];
		rs[i] = min;
    }
}






SEXP padjust( SEXP PVALUE )
{
	int j, i, N = LENGTH(PVALUE);
	int *lp = Calloc(N, int);
	int *ord = Calloc(N, int);
	int *reord = Calloc(N, int);
	double *ord2 = Calloc(N, double);
	double *ans = Calloc(N,double);
	SEXP s, pw;
	PROTECT(s = allocVector(REALSXP, LENGTH(PVALUE)));
	PROTECT(pw = allocVector(REALSXP, LENGTH(PVALUE)));
	for(i = N-1, j=0 ; i >= 0 ; i--, j++){
		ans[i] = -1*REAL(PVALUE)[i];
		REAL(s)[i] = NA_REAL;
		lp[j] = i+1;
	}
	for(i = 0 ; i < N ; i++){
		ord[i] = i;
		reord[i] = i;
	}
	if(N > 1){
		rsort_index(ans, ord, N);
	for(i = 0 ; i < N ; i++){
		ord2[i] = (double) ord[i];
	}
	rsort_index(ord2, reord, N);
	}
	for(i = 0 ; i < N ; i++){
		REAL(pw)[i] = ((double) N / (double) lp[i]) * REAL(PVALUE)[ord[i]];
	}
	cummin(pw,s);
	for(i = 0 ; i < N ; i++){
		if(1.0 < REAL(s)[reord[i]]){
			REAL(pw)[i] = 1.0;
		}else{
			REAL(pw)[i] = REAL(s)[reord[i]];
		}
	}
	Free(lp);Free(ord);Free(ord2);Free(reord);Free(ans);
	UNPROTECT(2);
	return(pw);
}


void padjust_C( double *PVALUE, int N )
{
	int j, i;
	int *lp = Calloc(N, int);
	int *ord = Calloc(N, int);
	int *reord = Calloc(N, int);
	double *ord2 = Calloc(N, double);
	double *ans = Calloc(N,double);
	double *s = Calloc(N, double);
	double *pw = Calloc(N, double);
	
	for(i = N-1, j=0 ; i >= 0 ; i--, j++){
		ans[i] = -1*PVALUE[i];
		s[i] = NA_REAL;
		lp[j] = i+1;
	}
	for(i = 0 ; i < N ; i++){
		ord[i] = i;
		reord[i] = i;
	}
	if(N > 1){
		rsort_index(ans, ord, N);
		for(i = 0 ; i < N ; i++){
			ord2[i] = (double) ord[i];
		}
		rsort_index(ord2, reord, N);
	}
	for(i = 0 ; i < N ; i++){
		pw[i] = ((double) N / (double) lp[i]) * PVALUE[ord[i]];
	}
	cummin_C(pw,s,N);
	for(i = 0 ; i < N ; i++){
		if(1.0 < s[reord[i]]){
			pw[i] = 1.0;
		}else{
			pw[i] = s[reord[i]];
		}
	}
	Free(lp);Free(ord);Free(ord2);Free(reord);Free(ans);
	for(i = 0 ; i < N ; i++){
		PVALUE[i] = pw[i];
	}
	Free(pw);Free(s);
}




SEXP split_summary( SEXP TBASE, SEXP tol )
{
    int i,j,k,l,N_splits = 0;
    double xx = 0.0;
    int N = LENGTH(TBASE);
	int *n_splits = Calloc(N, int);
    for(i=0; i < N; i++){
		/* the number of cuts per variable */
        n_splits[i] = LENGTH(VECTOR_ELT(VECTOR_ELT(TBASE,i),0));
        N_splits += n_splits[i];
    }

	int *Var_id = Calloc(N_splits, int);
	int *id_Var = Calloc(N_splits, int);
	double *X = Calloc(N_splits, double);

    k=0;
    for(j=0; j < N; j++){
        for(i=0; i < n_splits[j]; i++){
            X[k]= REAL(VECTOR_ELT(VECTOR_ELT(TBASE,j),0))[i];
            k++;
        }
    }
    l=0;
    for(i=0; i < N; i++){
        if(n_splits[i] != 0){
            for(j=1; j <= n_splits[i]; j++){
                Var_id[l]=j;
                id_Var[l]=i+1;
                l++;
            }
        }
    }
	Free(n_splits);
	l=0;
	xx = dmax(X,k);
	for(i=0; i < k; i++){
		if(X[i] >= xx*REAL(tol)[0]){
			Var_id[l]=Var_id[i];
			id_Var[l]=id_Var[i];
			X[l]=X[i];
			l++;
			}
		}
	k=l;
	int *i_in = Calloc(k, int);
	for(i=0; i < k; i++){
		i_in[i]=i;
	}
	rsort_with_index(X, i_in, k);
	SEXP X_in = PROTECT(allocVector(REALSXP,k));
	SEXP Var_id2 = PROTECT(allocVector(INTSXP,k));
	SEXP id_Var2 = PROTECT(allocVector(INTSXP,k));
	l=0;
	for(i=N_splits-1; i > (N_splits-1) - k; i--){
		INTEGER(Var_id2)[l]=Var_id[i_in[i]];
		INTEGER(id_Var2)[l]=id_Var[i_in[i]];
		REAL(X_in)[l]=X[i];
		l++;
	}
	Free(Var_id); Free(id_Var); Free(X); 
	Free(i_in);
	
	SEXP globD = PROTECT(allocVector(REALSXP,1));
    REAL(globD)[0] = REAL(VECTOR_ELT(VECTOR_ELT(TBASE,0),1))[0];
    SEXP which = PROTECT(allocVector(VECSXP,N));
    for(i=0; i < N; i++){
        SET_VECTOR_ELT(which, i, VECTOR_ELT(VECTOR_ELT(TBASE,i),2));
    }
    SEXP result_out = PROTECT(allocVector(VECSXP,5));
    SET_VECTOR_ELT(result_out, 0, X_in);
    SET_VECTOR_ELT(result_out, 1, globD);
    SET_VECTOR_ELT(result_out, 2, which);
    SET_VECTOR_ELT(result_out, 3, id_Var2);
    SET_VECTOR_ELT(result_out, 4, Var_id2);
    UNPROTECT(6);
    return(result_out);  
}



SEXP split_summary_padj( SEXP SPLITS, SEXP tol )
{
	int i, j, k, n, N_splits = 0;
	int N = LENGTH(SPLITS);
	int *n_splits = Calloc(N, int);
	/* in padj case - sort X with respect to statistic*/
	/* the number of cuts per variable */
	for(i=0; i < N; i++){
        n_splits[i] = LENGTH(VECTOR_ELT(VECTOR_ELT(SPLITS,i),0));
        N_splits += n_splits[i];
    }
	int *ind_split = Calloc(N_splits, int);
	int *ind_var = Calloc(N_splits, int);
	double *PVALUE = Calloc(N_splits, double);
	double *Stats = Calloc(N_splits, double);
	
    k=0;
    for(j=0; j < N; j++){
        for(i=0; i < n_splits[j]; i++){
            PVALUE[k] = REAL(VECTOR_ELT(VECTOR_ELT(SPLITS,j),0))[i];
            Stats[k] = REAL(VECTOR_ELT(VECTOR_ELT(SPLITS,j),3))[i];
			ind_split[k] = i+1;
			ind_var[k] = j+1;
            k++;
        }
    }
	Free(n_splits);
	int *ind_in = Calloc(N_splits, int);
	for(i=0; i < N_splits; i++){
		ind_in[i]=i;
	}
	rsort_xyz(PVALUE, Stats, ind_in, N_splits);
	n = 0;
	for(i=0; i < N_splits; i++){
		if(fabs(PVALUE[i]) <= FLT_EPSILON){
			Stats[i] = -1.0*Stats[i];
			n++;
		}
	}
	if(n > 1){
		rsort_with_index(Stats, ind_in, n);
	}
	padjust_C(PVALUE, k);
	
	Free(Stats);
	SEXP X = PROTECT(allocVector(REALSXP,k));
	SEXP IND_VAR = PROTECT(allocVector(INTSXP,k));
	SEXP IND_SPL = PROTECT(allocVector(INTSXP,k));
	j=0;
	for(i=0; i < k; i++){
		INTEGER(IND_VAR)[j] = ind_var[ind_in[i]];
		INTEGER(IND_SPL)[j] = ind_split[ind_in[i]];
		REAL(X)[j] = 1.0-PVALUE[i];
		j++;
	}
	Free(ind_var); Free(ind_split); Free(PVALUE);Free(ind_in);
	SEXP globD = PROTECT(allocVector(REALSXP,1));
    REAL(globD)[0] = REAL(VECTOR_ELT(VECTOR_ELT(SPLITS,0),1))[0];
    SEXP which = PROTECT(allocVector(VECSXP,N));
    for(i=0; i < N; i++){
        SET_VECTOR_ELT(which, i, VECTOR_ELT(VECTOR_ELT(SPLITS,i),2));
    }
    SEXP result_out = PROTECT(allocVector(VECSXP,5));
    SET_VECTOR_ELT(result_out, 0, X);
    SET_VECTOR_ELT(result_out, 1, globD);
    SET_VECTOR_ELT(result_out, 2, which);
    SET_VECTOR_ELT(result_out, 3, IND_VAR);
    SET_VECTOR_ELT(result_out, 4, IND_SPL);
    UNPROTECT(6);
    return(result_out); 
}





SEXP split_summary_dev( SEXP SPLITS, SEXP tol )
{
	int i, j, k, N_splits = 0;
	int N = LENGTH(SPLITS);
	int *n_splits = Calloc(N, int);
	
	/* the number of cuts per variable */
	for(i=0; i < N; i++){
        n_splits[i] = LENGTH(VECTOR_ELT(VECTOR_ELT(SPLITS,i),0));
        N_splits += n_splits[i];
    }
	int *ind_split = Calloc(N_splits, int);
	int *ind_var = Calloc(N_splits, int);
	double *DEV = Calloc(N_splits, double);
	
    k=0;
    for(j=0; j < N; j++){
        for(i=0; i < n_splits[j]; i++){
            DEV[k] = REAL(VECTOR_ELT(VECTOR_ELT(SPLITS,j),0))[i];
			ind_split[k] = i+1;
			ind_var[k] = j+1;
            k++;
        }
    }
	Free(n_splits);
	double max_DEV =0.0;
	max_DEV = dmax(DEV, N_splits);
	int *ind_in = Calloc(N_splits, int);
	for(i=0; i < N_splits; i++){
		ind_in[i]=i;
	}
	rsort_with_index(DEV, ind_in, N_splits);
	k = 0;
	for(i=0; i < N_splits; i++){
		if(DEV[i] >= max_DEV * REAL(tol)[0]){
			k++;
		}
	}
	SEXP X = PROTECT(allocVector(REALSXP,k));
	SEXP IND_VAR = PROTECT(allocVector(INTSXP,k));
	SEXP IND_SPL = PROTECT(allocVector(INTSXP,k));
	j=0;
	for(i=N_splits-1; i > (N_splits-1) - k; i--){
		INTEGER(IND_VAR)[j] = ind_var[ind_in[i]];
		INTEGER(IND_SPL)[j] = ind_split[ind_in[i]];
		REAL(X)[j] = DEV[i];
		j++;
	}
	/* result */
	SEXP globD = PROTECT(allocVector(REALSXP,1));
    REAL(globD)[0] = REAL(VECTOR_ELT(VECTOR_ELT(SPLITS,0),1))[0];
    SEXP which = PROTECT(allocVector(VECSXP,N));
    for(i=0; i < N; i++){
        SET_VECTOR_ELT(which, i, VECTOR_ELT(VECTOR_ELT(SPLITS,i),2));
    }
    SEXP result_out = PROTECT(allocVector(VECSXP,5));
    SET_VECTOR_ELT(result_out, 0, X);
    SET_VECTOR_ELT(result_out, 1, globD);
    SET_VECTOR_ELT(result_out, 2, which);
    SET_VECTOR_ELT(result_out, 3, IND_VAR);
    SET_VECTOR_ELT(result_out, 4, IND_SPL);
	Free(ind_var); Free(ind_split); Free(DEV);Free(ind_in);
    UNPROTECT(6);
    return(result_out); 
}



