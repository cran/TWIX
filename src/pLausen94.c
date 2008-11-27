#include <Rinternals.h>
#include <Rdefines.h>
#include <R.h>
#include <Rmath.h>

int d2i_round(double Zahl){
	return(Zahl<0?Zahl-.5:Zahl+.5);
}

int d2i_floor(double Zahl, int Stellen)
{
	double v[] = {1, 4, 1e2, 1e3, 1e4 };
	return(floor(Zahl * v[Stellen] + 0.5) / v[Stellen]);
}




SEXP pLausen94(SEXP q, SEXP N, SEXP minprop, SEXP maxprop, SEXP m)
{
	int i, k=0, N_m=0, start, end;
	double *NN = REAL(N);
	double *Q = REAL(q);
	SEXP result;
	
	if(isNull(m)){
		start = d2i_floor(REAL(N)[0]*REAL(minprop)[0],2);
		end = d2i_floor(REAL(N)[0]*REAL(maxprop)[0],2);
		PROTECT(m = allocVector(INTSXP,end-start+1));
		for (i = 0,k=start; i < end-start+1; i++){
			REAL(m)[i] = k;
			k++;
		}
		UNPROTECT(1);
	}
	N_m = LENGTH(m);
	double *m1, *m2, *T;
	m1 = Calloc(N_m, double);
	m2 = Calloc(N_m, double);
	T = Calloc(N_m-1, double);
	double *M = REAL(m);
	if(N_m < 2){
		m1[0] = M[0];
		m2[0] = M[0];
		N_m = 1;
	}
	else{
		for(i = 0; i < N_m-1; i++){
			m1[i] = M[i];
			m2[i] = M[i+1];
		}
	}
	/* compute t and D */

	double D = 0.0;
	for(i = 0; i < N_m-1; i++){
		T[i] = sqrt(1.0-(m1[i]*(NN[0]-m2[i]))/((NN[0]-m1[i])*m2[i]));
		D += (M_1_PI)*exp(-(pow(Q[0],2))/2)*(T[i] - ((pow(Q[0],2))/4 -1)*(pow(T[i],3))/6);
	}
	Free(m1);Free(m2);Free(T);

	PROTECT(result = allocVector(REALSXP,1));
	REAL(result)[0] = 1.0 - (pnorm(Q[0], 0.0, 1.0, 1, 0) - pnorm((-1.0)*Q[0], 0.0, 1.0, 1, 0)) + D;
	UNPROTECT(1);
	return(result);
}



double C_pLausen94( SEXP q, double NN, double *M, int N_m)
{
	int i;
	double *Q = REAL(q);
	
	double *m1, *m2, *T;
	m1 = Calloc(N_m, double);
	m2 = Calloc(N_m, double);
	T = Calloc(N_m-1, double);
	if(N_m < 2){
		m1[0] = M[0];
		m2[0] = M[0];
		N_m = 1;
	}
	else{
		for(i = 0; i < N_m-1; i++){
			m1[i] = M[i];
			m2[i] = M[i+1];
		}
	}
	/* compute t and D */
	double D = 0.0;
	for(i = 0; i < N_m-1; i++){
		T[i] = sqrt(1.0-(m1[i]*(NN-m2[i]))/((NN-m1[i])*m2[i]));
		D += (M_1_PI)*exp(-(pow(Q[0],2))/2)*(T[i] - ((pow(Q[0],2))/4 -1)*(pow(T[i],3))/6);
	}
	Free(m1);Free(m2);Free(T);
	return(1.0 - (pnorm(Q[0], 0.0, 1.0, 1, 0) - pnorm((-1.0)*Q[0], 0.0, 1.0, 1, 0)) + D);
}







