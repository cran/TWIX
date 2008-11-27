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


static double Devleaf2(double *x, int n){

	double result=0.0;
	double count,j,TD=0.0;
	int max_lev,i=0;
	
	max_lev=dmax(x,n);
	for(j=0.0; j < max_lev; j++){
		count=0.0;
		for(i=0; i < n; i++){
			if(x[i] == j+1.0){
				count=count+1.0;
			}
		}
		if(count/n != 0){
			TD = TD + count*log(count/n);
		}
	}
	result= -TD;
	return(result);
}


SEXP Dev_oob(SEXP obj, SEXP data, SEXP rsp){

	SEXP dims,Imp_Dev,levels,data_cat;
	int k=0,count=0,i=0,j=0;
	double GDEV=0.0;
	double *v_rsp=REAL(rsp);
	int N = LENGTH(rsp);
	int sp_vek[N];
	
	dims=getAttrib(data, R_DimSymbol);
	GDEV = Devleaf2(v_rsp,LENGTH(rsp));
	for(i=0; i < N; i++){
		sp_vek[i]=0;
	}
	if(!isVectorList(data)){
		SEXP L_rsp,R_rsp;
		PROTECT(Imp_Dev = allocVector(REALSXP,1));
		if(LENGTH(VECTOR_ELT(obj, 0)) > 0 && VECTOR_ELT(obj, 0) != 0){
			if(inherits(data,"factor")){
				levels = getAttrib(data, R_LevelsSymbol);
				data_cat=coerceVector(data,INTSXP);
				for(j=0,k=0; j < LENGTH(levels); j++,k+=2){
					if(INTEGER(VECTOR_ELT(obj, 2))[k] == 1){
						for(i=0; i < N; i++){
							if(INTEGER(data_cat)[i] == j+1){
								sp_vek[i] = 1;
								count++;
							}
						}
					}
				}
			}
			else{
				for(i=0; i < N; i++){
					if(REAL(data)[i] < REAL(VECTOR_ELT(obj, 2))[0]){
						sp_vek[i] = 1;
						count++;
					}
				}
			}
			if(count > 0){
				int r,l;				
				PROTECT(L_rsp = allocVector(REALSXP,count));
				PROTECT(R_rsp = allocVector(REALSXP,N-count));

 				for(i=0,r=0,l=0; i < N; i++){
					if(sp_vek[i] == 1){
						REAL(L_rsp)[l]=REAL(rsp)[i];
						l++;
					}
					else{
						REAL(R_rsp)[r]=REAL(rsp)[i];
						r++;
					}
				}
			}
			else{
				PROTECT(L_rsp = allocVector(REALSXP,1));
				PROTECT(R_rsp = allocVector(REALSXP,1));
				REAL(L_rsp)[0]=0;
				REAL(R_rsp)[0]=0;
			}
			double *lv_rsp=REAL(L_rsp);
			double *rv_rsp=REAL(R_rsp);
			REAL(Imp_Dev)[0]=GDEV-Devleaf2(lv_rsp,LENGTH(L_rsp))-Devleaf2(rv_rsp,LENGTH(R_rsp));
			UNPROTECT(2);
		}
		else{
			REAL(Imp_Dev)[0]=0;
		}
	}
	else{
		SEXP list_obj,data_obj;
		int NL,h;
		NL=LENGTH(obj);
		PROTECT(Imp_Dev = allocVector(REALSXP,NL));
		for(h=0; h < NL; h++){
			for(i=0; i < N; i++){
					sp_vek[i]=0;
			}
			i=0,count=0;
			list_obj=VECTOR_ELT(obj, h);
			data_obj=VECTOR_ELT(data, h);
			if(LENGTH(VECTOR_ELT(list_obj, 0)) > 0 && VECTOR_ELT(list_obj, 0) != 0){
				if(inherits(data_obj,"factor")){
					levels = getAttrib(data_obj, R_LevelsSymbol);
					data_obj=coerceVector(data_obj,INTSXP);
					if(LENGTH(VECTOR_ELT(list_obj, 2)) > 2){
						for(j=0,k=0; j < LENGTH(levels); j++,k+=2){
							if(INTEGER(VECTOR_ELT(list_obj, 2))[k] == 1){
								for(i=0; i < N; i++){
									if(INTEGER(data_obj)[i] == j+1){
										sp_vek[i] = 1;
										count++;
									}
								}
							}
						}
					}
					else{
						for(j=0,k=0; j < LENGTH(levels); j++,k+=1){
							if(INTEGER(VECTOR_ELT(list_obj, 2))[k] == 1){
								for(i=0; i < N; i++){
									if(INTEGER(data_obj)[i] == j+1){
										sp_vek[i] = 1;
										count++;
									}
								}
							}
						}
					}
				}
				else{
					data_obj=coerceVector(data_obj,REALSXP);
					for(i=0; i < N; i++){
						if(REAL(data_obj)[i] < REAL(VECTOR_ELT(list_obj, 2))[0]){
							sp_vek[i] = 1;
							count++;
						}
					}
				}
				double LL_rsp[count],RR_rsp[N-count];
				if(count > 0){
					int r,l;				
					// PROTECT(L_rsp = allocVector(REALSXP,count));
					// PROTECT(R_rsp = allocVector(REALSXP,N-count));
					
//					L_rsp = allocVector(REALSXP,count);
//					R_rsp = allocVector(REALSXP,N-count);

					for(i=0,r=0,l=0; i < N; i++){
						if(sp_vek[i] == 1){
//							REAL(L_rsp)[l]=REAL(rsp)[i];
							LL_rsp[l]=REAL(rsp)[i];
							l++;
						}
						else{
//							REAL(R_rsp)[r]=REAL(rsp)[i];
							RR_rsp[r]=REAL(rsp)[i];
							r++;
						}
					}
				}
//				double *lv_rsp=REAL(L_rsp);
//				double *rv_rsp=REAL(R_rsp);
//				REAL(Imp_Dev)[h]=GDEV-Devleaf2(lv_rsp,LENGTH(L_rsp))-Devleaf2(rv_rsp,LENGTH(R_rsp));
				REAL(Imp_Dev)[h]=GDEV-Devleaf2(LL_rsp,count)-Devleaf2(RR_rsp,N-count);
				// UNPROTECT(2);
			}
			else{
				REAL(Imp_Dev)[h]=0;
			}
		}
	}
	UNPROTECT(1);
	return(Imp_Dev);
}

