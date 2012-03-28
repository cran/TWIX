
#include "utils.h"



/******************************************
 asymptotic 1 - P-value proposed by Lausen 1994.
 *\param Q quantile
 *\param N a number of observations
 *\param m an integer vector containing the sample sizes of the first group
 *\param N_m length of the M vector
 */


double C_pLausen94(const double Q, double N, const double *m, int N_m)
{
	int i=0;
	double *m1 = Calloc(N_m, double);
	double *m2 = Calloc(N_m, double);
	double *T = Calloc(N_m-1, double);
	if(N_m < 2){
		m1[0] = m[0];
		m2[0] = m[0];
		N_m = 1;
	}
	else{
		for(i = 0; i < N_m-1; i++){
			m1[i] = m[i];
			m2[i] = m[i+1];
		}
	}
	/* compute t and D */
	double D = 0.0, pval = 0.0;

	for(i = 0; i < N_m-1; i++){
		T[i] = sqrt(1.0-(m1[i]*(N-m2[i]))/((N-m1[i])*m2[i]));
		D += (M_1_PI)*exp(-(pow(Q,2))/2)*(T[i] - ((pow(Q,2))/4 -1)*(pow(T[i],3))/6);
	}
	Free(m1);Free(m2);Free(T);
	pval = 2.0 * (1.0 - pnorm(Q, 0.0, 1.0, 1, 0)) + D;
		
	if(pval > 1.0){
		pval = 1.0;
	}
	if(pval <= FLT_EPSILON){
		pval = 0.0;
	}
	//*return(1.0 - pval);
	return(pval);
}



/******************************************
 asymptotic 1 - P-value proposed by Lausen 1994.
 *\param Q a numeric vector of quantiles
 *\param N a number of observations
 *\param m an integer vector containing the sample sizes of the first group
 *\param N_m length of the M vector
 *\param pval a vector of 1 - p-values
 */


void C_pLausen94_all(const double *Q, double N, const double *m, int N_m, double *pval)
{
	int i,j;
	double *m1 = Calloc(N_m, double);
	double *m2 = Calloc(N_m, double);
	double *T = Calloc(N_m-1, double);
	if(N_m < 2){
		m1[0] = m[0];
		m2[0] = m[0];
		N_m = 1;
	}
	else{
		for(i = 0; i < N_m-1; i++){
			m1[i] = m[i];
			m2[i] = m[i+1];
		}
	}
	/* compute t and D */
	for(j = 0; j < N_m; j++){
		pval[j] = 0.0;
		double D = 0.0;
		for(i = 0; i < N_m-1; i++){
			T[i] = sqrt(1.0-(m1[i]*(N-m2[i]))/((N-m1[i])*m2[i]));
			D += (M_1_PI)*exp(-(pow(Q[j],2))/2)*(T[i] - ((pow(Q[j],2))/4 -1)*(pow(T[i],3))/6);
		}
		pval[j] = 2.0 * (1.0 - pnorm(Q[j], 0.0, 1.0, 1, 0)) + D;
		if(pval[j] > 1.0){
			pval[j] = 1.0;
		}
		if(pval[j] <= FLT_EPSILON){
			pval[j] = 0.0;
		}
		//*pval[j] = 1.0 - pval[j];
	}
	if(N_m - 1 < 1){
		pval[0] = 2.0 * (1.0 - pnorm(Q[0], 0.0, 1.0, 1, 0));
		if(pval[0] > 1.0){
			pval[0] = 1.0;
		}
		if(pval[0] <= FLT_EPSILON){
			pval[0] = 0.0;
		}
		//*pval[0] = 1.0 - pval[0];
	}
	Free(m1);Free(m2);Free(T);
}



/******************************************
 asymptotic 1 - P-value proposed by Lausen 1992.
 *\param Q  quantile
 *\param minp a minimum proportion of the observations
 *\param maxp a maximum proportion of the observations
*/


double C_pLausen92(const double Q, double minp, double maxp)
{
	double pval = 0.0, DQ = 0.0;
	if(Q <= 1.0){
		return(1.0);
	}
	DQ = dnorm4(Q, 0.0, 1.0, 0);
	pval = 4.0 * DQ / Q + DQ*(Q - 1.0/Q)*log((maxp*(1.0 - minp))/((1.0 - maxp)*minp));
//*	return(1.0 - pval);
	return(pval);
}



/******************************************
 asymptotic 1 - P-value proposed by Lausen 1992.
 *\param Q  a numeric vector of quantiles
 *\param N length of vector Q
 *\param minp a minimum proportion of the observations
 *\param maxp a maximum proportion of the observations
 *\param pval a vector of 1 - p-values
 */


void C_pLausen92_all(const double *Q, int N, double minp, double maxp, double *pval)
{
	double DQ;
	int i;
	for(i=0; i<N; i++){
		pval[i] = 1.0;
		DQ = 0.0;
		if(Q[i] > 1.0){
			DQ = dnorm4(Q[i], 0.0, 1.0, 0);
//*			pval[i] = 1.0 - (4.0 * DQ / Q[i] + DQ*(Q[i] - 1./Q[i])*log((maxp*(1. - minp))/((1.-maxp)*minp)));
			pval[i] = (4.0 * DQ / Q[i] + DQ*(Q[i] - 1.0/Q[i])*log((maxp*(1.0 - minp))/((1.0-maxp)*minp)));
		}
	}
}




