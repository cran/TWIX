

#include "utils.h"



int d2i_round(double Zahl)
{
	return(Zahl<0?Zahl-.5:Zahl+.5);
}



int d2i_floor(double Zahl, int Stellen)
{
	double v[] = {1, 4, 1e2, 1e3, 1e4 };
	return(floor(Zahl * v[Stellen] + 0.5) / v[Stellen]);
}



double clogn(double x){
	if (x <= FLT_EPSILON) 
		return(0.0); 
	else 
		return(x*log(x));
}



int graycode(int *indy,  int x){
	int i;
    int maxc=x;

    for( i=0; i < (maxc-1); i++) {
        if(indy[i] == 1) {
			indy[i]=2;
			return(i);
        }
		else if(indy[i] == 2) indy[i]=1;
    }
    return(maxc);
}



double dmax(double *X, int n)
{
	int i=0;
	double max=0.0;
	for (i=0; i < n; i++){
		if (X[i] > max) {
			max=X[i];
		}
	}
	return(max);
}




double SUMV_D(double *x, int N)
{
	int i;
	LDOUBLE sum=0.0;
    for (i = 0; i < N; i++){
		sum += x[i];
	}
    return(sum);
}




int SUMV_I(int *x, int N)
{
    int i;
	int sum=0;
    for (i = 0; i < N; i++){
		sum += x[i];
	}
    return(sum);
}





void cumsum(const double *x, int nx, double *ans)
{
    int i;
    LDOUBLE sum = 0.;
    for (i = 0 ; i < nx ; i++) {
		sum += x[i];
		ans[i] = sum;
    }
}




static int rcmp_TW(double x, double y, Rboolean nalast)
{
    int nax = ISNAN(x), nay = ISNAN(y);
    if (nax && nay)	return 0;
    if (nax)		return nalast?1:-1;
    if (nay)		return nalast?-1:1;
    if (x < y)		return -1;
    if (x > y)		return 1;
    return 0;
}



void rsort_index(double *x, int *indx, int n)
{
    double v;
    int i, j, h, iv;

    for (h = 1; h <= n / 9; h = 3 * h + 1);
    for (; h > 0; h /= 3)
	for (i = h; i < n; i++) {
	    v = x[i]; iv = indx[i];
	    j = i;
	    while (j >= h && rcmp_TW(x[j - h], v, TRUE) > 0)
		 { x[j] = x[j - h]; indx[j] = indx[j-h]; j -= h; }
	    x[j] = v; indx[j] = iv;
	}
}




void c_rsort(double *x, int n)
{
    double v;
    int i, j, h;

    for (h = 1; h <= n / 9; h = 3 * h + 1);
    for (; h > 0; h /= 3)
	for (i = h; i < n; i++) {
	    v = x[i];
	    j = i;
	    while (j >= h && rcmp_TW(x[j - h], v, TRUE) > 0)
		 { x[j] = x[j - h]; j -= h; }
	    x[j] = v;
	}
}



void rsort_with_x(double *x, double *indx, int n)
{
    double v, iv;
    int i, j, h;
	
    for (h = 1; h <= n / 9; h = 3 * h + 1);
    for (; h > 0; h /= 3)
	for (i = h; i < n; i++) {
	    v = x[i]; iv = indx[i];
	    j = i;
	    while (j >= h && rcmp_TW(x[j - h], v, TRUE) > 0)
		 { x[j] = x[j - h]; indx[j] = indx[j-h]; j -= h; }
		x[j] = v; indx[j] = iv;
	}
}


void rsort_xyz(double *x, double *y, int *indx, int n)
{
    double v, iv, vi;
    int i, j, h;
	
    for (h = 1; h <= n / 9; h = 3 * h + 1);
    for (; h > 0; h /= 3)
		for (i = h; i < n; i++) {
			v = x[i]; iv = indx[i]; vi = y[i];
			j = i;
			while (j >= h && rcmp_TW(x[j - h], v, TRUE) > 0)
			{ x[j] = x[j - h]; indx[j] = indx[j-h]; y[j] = y[j-h]; j -= h; }
			x[j] = v; indx[j] = iv; y[j] = vi;
		}
}


int nrow(SEXP x) {
    return(INTEGER(getAttrib(x, R_DimSymbol))[0]);
}


int ncol(SEXP x) {
    return(INTEGER(getAttrib(x, R_DimSymbol))[1]);
}



// X is a data.frame

SEXP mysplit(SEXP X){

	int i, j;
	int n = LENGTH(VECTOR_ELT(X, 0));
	int m = LENGTH(X);
	
	SEXP ans = PROTECT(allocVector(VECSXP, n));
	for(i = 0; i < n; i++) {
		SEXP C_DATA = PROTECT(allocVector(VECSXP,m));
		for(j = 0; j < m; j++){
			int mode = TYPEOF(VECTOR_ELT(X, j));
			SEXP C_NRow = PROTECT(allocVector(mode,1));
			switch (mode){
				case LGLSXP:
				case INTSXP:
						INTEGER(C_NRow)[0] = INTEGER(VECTOR_ELT(X, j))[i];
					break;
				case REALSXP:
						REAL(C_NRow)[0] = REAL(VECTOR_ELT(X, j))[i];
					break;
				case STRSXP:
						SET_STRING_ELT(C_NRow,0,STRING_ELT(VECTOR_ELT(X, j),i));
					break;
			}
			if(isFactor(VECTOR_ELT(X, j))){
				DUPLICATE_ATTRIB(C_NRow,VECTOR_ELT(X, j));
			}
			SET_VECTOR_ELT(C_DATA, j, C_NRow);
			UNPROTECT(1);
		}
/*		
	dim-names if X is a matrix
		
		PROTECT(xdims = allocVector(INTSXP, 2));
		INTEGER(xdims)[0] = m;
		INTEGER(xdims)[1] = 1;
		setAttrib(DATA,R_DimSymbol,xdims);
		UNPROTECT(1);
*/		
		SEXP row_names = PROTECT(allocVector(INTSXP, 2));
		INTEGER(row_names)[0] = NA_INTEGER;
		INTEGER(row_names)[1] = 1;
		setAttrib(C_DATA, R_RowNamesSymbol, row_names);
		UNPROTECT(1);
		
		namesgets(C_DATA,getAttrib(X,R_NamesSymbol));
		setAttrib(C_DATA, R_ClassSymbol, mkString("data.frame"));
		SET_VECTOR_ELT(ans, i, C_DATA);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return ans;
}


/******************************************
* search all of the local maxima of the x
* 
* @param:x is vector 
* @param:N is length of the x
* @param:indx index of the new sequence
*/

void loc_maxima( double *x, int N, int *indx )
{
	double x_max = 0.0;
	int i, count = 0;
	
	for(i=0; i < N-2; i++){
		if(x[i+1] > x[i] && x[i+1] > x[i+2]){
			indx[i+1] = 1;
			count++;
		}
	}
	x_max = dmax(x,N);
	if(fabs(x[0] - x_max) <= FLT_EPSILON){
		indx[0]=1;
		count++;
	}
	if(fabs(x[N-1] - x_max) <= FLT_EPSILON){
		indx[N-1]=1;
		count++;
	}			
}


 
/******************************************
* sequence of the x
* 
* @param:x is vector 
* @param:N is length of the x
* @param:step increment of the sequence
* @param:indx index of the new sequence
*/

void dev_seq( int N, int step, int *indx )
{
	int from = 0;
	int to = N-1;
	int i, j=0, del=0; 
	double m=0.0;
	
	del = to - from;
	m =  del / step;
	j = (int) m;

	for(i=0; i < j+1; i++){
		indx[i*(step)] = 1;
	}
}





/******************************************
* combination methods
*
* @param : pred  list of predictions(matrix) 
* @param : certanty a vector of weights
*/

SEXP bayes_comb( SEXP pred, SEXP certanty)
{
	// number of prediction matrices
	int N = LENGTH(pred);
	int nobs = LENGTH(VECTOR_ELT(pred,0));
		
	int i, j;
	for(j=0; j < nobs; j++){
		REAL(VECTOR_ELT(pred,0))[j] = REAL(certanty)[0] * REAL(VECTOR_ELT(pred,0))[j];
	}
	
	for(j=0; j < nobs; j++){
		for(i=1; i < N; i++){
			REAL(VECTOR_ELT(pred,0))[j] += REAL(certanty)[i] * REAL(VECTOR_ELT(pred,i))[j];
		}
	}
	return(pred);
}



