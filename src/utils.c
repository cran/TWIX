#include <Rinternals.h>
#include <Rdefines.h>
#include <R.h>
#include <Rmath.h>



double dmax(double *X,int n)
{
	int i=0;
	double max=0.0;
	for (i=0;i < n; i++) {
		if (X[i] > max) {
			max=X[i];
		}
	}
	return(max);
}

double SUMV_D(double *x, int N)
{
	int i;
	double sum=0.0;
    for (i = 0; i < N; i++){
		sum += x[i];
	}
    return sum;
}

int SUMV_I(int *x, int N)
{
    int i;
	int sum=0;
    for (i = 0; i < N; i++){
		sum += x[i];
	}
    return sum;
}


void cumsum(double *x, int nx, double *ans)
{
    int i;
    double sum = 0.0;
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
