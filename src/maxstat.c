#include <Rinternals.h>
#include <Rdefines.h>
#include <R.h>
#include <Rmath.h>

#include "pLausen94.h"
#include "utils.h"


SEXP maxstat(SEXP Y, SEXP X, SEXP minprop, SEXP maxprop)
{
	SEXP result, STATISTIC, ESTIMATOR, STATS, CUTS;
	SEXP names_result, PVAL, ties, X_SORT;
	int i, k=0;
	double SS = 0.0;
	int N = LENGTH(Y), NM = 0;
	double N_d = N;
	
	PROTECT(X_SORT = allocVector(REALSXP,N));

	int *m_i;
	double *YY, *SY, *m_d;
	m_i = Calloc(N, int);
	m_d = Calloc(N, double);
	YY = Calloc(N, double);
	SY = Calloc(N, double);

	for (i = 0; i < N; i++){
		REAL(X_SORT)[i] = REAL(X)[i];
		YY[i] = REAL(Y)[i];
	}

	double *XX = REAL(X_SORT);
	rsort_with_x(XX,YY,N);

	ties = duplicated(X_SORT,0);
	UNPROTECT(1);

	for(i=d2i_floor(N_d*REAL(minprop)[0],2),k=0; i < d2i_floor(N_d*REAL(maxprop)[0],2)+1; i++){
		if(!LOGICAL(ties)[i]){
			m_i[k]=i;
			m_d[k]=i;
			k++;
		}
	}
	NM=k;
	if(NM > 0 && N > 1){
		double SSY, *E, *V, *Test;
		E = Calloc(NM, double);
		V = Calloc(NM, double);
		Test = Calloc(NM, double);
	
		SS = SUMV_D(YY,N);
		/* sum(y^2) */
		for (i = 0; i < N; i++){
			SY[i] = pow(YY[i],2);
		}
		SSY = SUMV_D(SY,N);
		/* cumsum(y) */
		cumsum(YY,N,SY);
		/* Test */
		for (i = 0; i < NM; i++){
			E[i] = (m_d[i])/N_d*SS;
			V[i] = m_d[i]*(N_d-m_d[i])/(pow(N_d,2)*(N_d-1.0))*(N_d*SSY-pow(SS,2));
			Test[i] = fabs((SY[m_i[i]-1] - E[i])/sqrt(V[i]));
		}
		Free(E);Free(V);Free(SY);Free(YY);
		PROTECT(STATISTIC = allocVector(REALSXP,1));
		REAL(STATISTIC)[0] = dmax(Test,NM);
		PROTECT(STATS = allocVector(REALSXP,NM));
		PROTECT(CUTS = allocVector(REALSXP,NM));
		PROTECT(ESTIMATOR = allocVector(REALSXP,1));
		REAL(ESTIMATOR)[0] = NA_REAL;
		SEXP which;
		PROTECT(which = allocVector(REALSXP,1));
		REAL(which)[0] = NA_REAL;
		/* ESTIMATOR */
		int update = 1;	
		for (i = 0; i < NM; i++){
			REAL(STATS)[i] = Test[i];
			REAL(CUTS)[i] = XX[m_i[i]-1];
			if(Test[i] == REAL(STATISTIC)[0]){
				if(update){
					REAL(ESTIMATOR)[0] = XX[m_i[i]-1];
					REAL(which)[0] = (XX[m_i[i]-1]+XX[m_i[i]])/2;
					update = 0;
				}
			}
		}
		Free(Test);Free(m_i);
		/* OUTPUT */
		PROTECT(names_result = allocVector(STRSXP, 6));
		SET_STRING_ELT(names_result, 0, mkChar("statistic"));
		SET_STRING_ELT(names_result, 1, mkChar("p.value"));
		SET_STRING_ELT(names_result, 2, mkChar("estimate"));
		SET_STRING_ELT(names_result, 3, mkChar("stats"));
		SET_STRING_ELT(names_result, 4, mkChar("cuts"));
		SET_STRING_ELT(names_result, 5, mkChar("which"));
		
		PROTECT(result = allocVector(VECSXP, 6));
		PROTECT(PVAL = allocVector(REALSXP,1));
		REAL(PVAL)[0] = C_pLausen94(STATISTIC,N_d,m_d,NM);
		Free(m_d);
		if(REAL(PVAL)[0] > 1)
			REAL(PVAL)[0] = 1.0;
	
		SET_VECTOR_ELT(result, 0, STATISTIC);
		SET_VECTOR_ELT(result, 1, PVAL);
		SET_VECTOR_ELT(result, 2, ESTIMATOR);
		SET_VECTOR_ELT(result, 3, STATS);
		SET_VECTOR_ELT(result, 4, CUTS);
		SET_VECTOR_ELT(result, 5, which);
	
		setAttrib(result, R_NamesSymbol, names_result);
		UNPROTECT(8);
		return(result);
	}
	else{
		PROTECT(names_result = allocVector(STRSXP, 5));
		SET_STRING_ELT(names_result, 0, mkChar("statistic"));
		SET_STRING_ELT(names_result, 1, mkChar("p.value"));
		SET_STRING_ELT(names_result, 2, mkChar("estimate"));
		SET_STRING_ELT(names_result, 3, mkChar("stats"));
		SET_STRING_ELT(names_result, 4, mkChar("cuts"));
		
		PROTECT(result = allocVector(VECSXP, 5));
		PROTECT(PVAL = allocVector(REALSXP,1));
		REAL(PVAL)[0] = NA_REAL;
		
		SET_VECTOR_ELT(result, 0, PVAL);
		SET_VECTOR_ELT(result, 1, PVAL);
		SET_VECTOR_ELT(result, 2, PVAL);
		SET_VECTOR_ELT(result, 3, PVAL);
		SET_VECTOR_ELT(result, 4, PVAL);
	
		setAttrib(result, R_NamesSymbol, names_result);
		UNPROTECT(3);
		return(result);
	}
}


