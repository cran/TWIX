
#include "utils.h"
#include "pLausen94.h"



SEXP maxstat(SEXP X, SEXP Y, SEXP minprop, SEXP maxprop, SEXP test, SEXP minbuck)
{
	SEXP result;
	SEXP names_result, PVAL, globD;
	int i, k=0;
	int N = LENGTH(Y), NM = 0;
	double N_d = N;
	
	int *m_i = Calloc(N, int);
	double *m_d = Calloc(N, double);
	double *YY = Calloc(N, double);
	double *X_order = Calloc(N, double);
	
	for (i = 0; i < N; i++){
		X_order[i] = REAL(X)[i];
		YY[i] = INTEGER(Y)[i];
	}
	rsort_with_x(X_order, YY, N);
	int start_inx = d2i_floor(N_d*REAL(minprop)[0],2);
	int end_inx = d2i_floor(N_d*REAL(maxprop)[0],2);
	
	if(start_inx < INTEGER(minbuck)[0]-1){
		start_inx = INTEGER(minbuck)[0]-1;
	}
	if(end_inx > N - INTEGER(minbuck)[0]){
		end_inx = N - INTEGER(minbuck)[0];
	}
	int *ties = Calloc(N, int);
	for(i=1; i < N; i++){
		if(fabs(X_order[i]-X_order[i-1]) <= FLT_EPSILON){
			ties[i] = 1;
		}
	}
	double Xtmp_last = X_order[0];
	for(i=start_inx; i <= end_inx; i++){
		if(ties[i] != 1){
			m_i[k]=i;
			m_d[k]=i;
			k++;
		}
	}
	Free(ties);
	NM=k;
	if(NM > 0 && N > 1){
		double *E, *V, *Test, *SY;
		double SSY= 0.0, SS = 0.0;

		SY = Calloc(N, double);
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
		
		/* Teststat */
		for (i = 0; i < NM; i++){
			E[i] = (m_d[i])/N_d*SS;
			V[i] = m_d[i]*(N_d-m_d[i])/(pow(N_d,2)*(N_d-1.0))*(N_d*SSY-pow(SS,2));
			if(V[i] > FLT_EPSILON){
				Test[i] = fabs((SY[m_i[i]-1] - E[i])/sqrt(V[i]));
			}else{
				Test[i] = 0.0;
			}
		}
		Free(E);Free(V);Free(SY);Free(YY);
		
		SEXP STATISTIC, ESTIMATOR, STATS, CUTS;
		
		PROTECT(STATISTIC = allocVector(REALSXP,1));
		REAL(STATISTIC)[0] = dmax(Test,NM);
		
		PROTECT(STATS = allocVector(REALSXP,NM));
		PROTECT(CUTS = allocVector(REALSXP,NM));
		PROTECT(ESTIMATOR = allocVector(REALSXP,1));
		REAL(ESTIMATOR)[0] = 0.0;
		SEXP which = PROTECT(allocVector(REALSXP,1));
		REAL(which)[0] = 0.0;
		/* ESTIMATOR */
		int update = 1;	
		for (i = 0; i < NM; i++){
			REAL(STATS)[i] = Test[i];
			REAL(CUTS)[i] = X_order[m_i[i]-1];
//*			Rprintf("\n i=%d, k=%d, Test[i]=%f, CUTS = %f, m_i[i] = %d, X_order[m_i[i]-1] = %f, X_order[m_i[i]] = %f \n",i,k,Test[i],REAL(CUTS)[i], m_i[i], X_order[m_i[i]-1],X_order[m_i[i]]);
			if(fabs(Test[i] - REAL(STATISTIC)[0]) <= FLT_EPSILON){
				if(update){
					REAL(ESTIMATOR)[0] = X_order[m_i[i]-1];
					REAL(which)[0] = (X_order[m_i[i]-1]+X_order[m_i[i]])/2.0;
					update = 0;
				}
			}
		}
		Free(m_i); Free(X_order);
		/* OUTPUT */
		if(LOGICAL(test)[0]){
			PROTECT(names_result = allocVector(STRSXP, 6));
			SET_STRING_ELT(names_result, 0, mkChar("statistic"));
			SET_STRING_ELT(names_result, 1, mkChar("p.value"));
			SET_STRING_ELT(names_result, 2, mkChar("estimate"));
			SET_STRING_ELT(names_result, 3, mkChar("stats"));
			SET_STRING_ELT(names_result, 4, mkChar("cuts"));
			SET_STRING_ELT(names_result, 5, mkChar("which"));
		
			PROTECT(result = allocVector(VECSXP, 6));
			if(NM < 50){
				if(NM==1){
					PROTECT(PVAL = allocVector(REALSXP,1));
				}else{
					PROTECT(PVAL = allocVector(REALSXP,NM));
				}
				C_pLausen94_all(Test, N_d, m_d, NM, REAL(PVAL));
			}
			else{
				PROTECT(PVAL = allocVector(REALSXP,NM));
				C_pLausen92_all(Test, NM, REAL(minprop)[0], REAL(maxprop)[0], REAL(PVAL));
			}
			Free(m_d); Free(Test);
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
			SEXP globD;
			PROTECT(names_result = allocVector(STRSXP, 4));
			SET_STRING_ELT(names_result, 0, mkChar("dev"));
			SET_STRING_ELT(names_result, 1, mkChar("globD"));
			SET_STRING_ELT(names_result, 2, mkChar("which"));
			SET_STRING_ELT(names_result, 3, mkChar("stats"));
			PROTECT(result = allocVector(VECSXP, 4));
			PROTECT(PVAL = allocVector(REALSXP,1));
			PROTECT(globD = allocVector(REALSXP,1));
			if(NM < 50){
				REAL(PVAL)[0] = C_pLausen94(REAL(STATISTIC)[0],N_d,m_d,NM);
			}
			else{ 
				REAL(PVAL)[0] = C_pLausen92(REAL(STATISTIC)[0],REAL(minprop)[0],REAL(maxprop)[0]);
			}
			Free(m_d);Free(Test);
			
			SET_VECTOR_ELT(result, 0, PVAL);
			REAL(globD)[0] = 1.0;
			SET_VECTOR_ELT(result, 1, globD);
			SET_VECTOR_ELT(result, 2, which);
			SET_VECTOR_ELT(result, 3, STATISTIC);
			setAttrib(result, R_NamesSymbol, names_result);
			UNPROTECT(9);
			return(result);
		}
	}
	else{
		PROTECT(names_result = allocVector(STRSXP, 4));
		SET_STRING_ELT(names_result, 0, mkChar("dev"));
		SET_STRING_ELT(names_result, 1, mkChar("globD"));
		SET_STRING_ELT(names_result, 2, mkChar("which"));
		SET_STRING_ELT(names_result, 3, mkChar("stats"));
		
		PROTECT(result = allocVector(VECSXP, 4));
		PROTECT(PVAL = allocVector(REALSXP,1));
		REAL(PVAL)[0] = 1.0;
		SEXP STATISTIC;
		PROTECT(STATISTIC = allocVector(REALSXP,1));
		REAL(STATISTIC)[0] = 0.0;
		SET_VECTOR_ELT(result, 0, PVAL);
		PROTECT(globD = allocVector(REALSXP,1));
		REAL(globD)[0] = 1.0;
		SET_VECTOR_ELT(result, 1, globD);
		SET_VECTOR_ELT(result, 2, PVAL);
		SET_VECTOR_ELT(result, 3, STATISTIC);
		setAttrib(result, R_NamesSymbol, names_result);
		UNPROTECT(5);
		Free(m_i);Free(m_d);Free(YY);Free(X_order);
		return(result);
	}
}


