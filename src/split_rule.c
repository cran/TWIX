#include <Rinternals.h>
#include <Rdefines.h>
#include <R.h>
#include <Rmath.h>

SEXP split_rule( SEXP DEV, SEXP GDEV, SEXP WHICH, SEXP MINBUCK,
                SEXP MINDEV, SEXP DATA, SEXP TDATA )
{
    SEXP result_out,Ldata,Ltdata,ans;
    double dev = REAL(DEV)[0];
    int minbuck = INTEGER(MINBUCK)[0];
    double mindev = REAL(MINDEV)[0];
    int NWh = LENGTH(WHICH);
    int ND = LENGTH(DATA);
    int NTD = LENGTH(TDATA);
    int i,j,k,Lobs=0,Robs=0,sw=0,sob=0;
    
    PROTECT(Ldata = allocVector(INTSXP,ND));
    PROTECT(Ltdata = allocVector(INTSXP,NTD));
    for(i=0; i < ND; i++){
        INTEGER(Ldata)[i]=0;
        }
    for(i=0; i < NTD; i++){
        INTEGER(Ltdata)[i]=0;
        }
    if( LENGTH(DEV) != 0) {
        if(NWh == 1) {
            if(dev != 0.0 ){
                sw = (dev != 0.0) && (dev/REAL(GDEV)[0] > mindev);
            }
            else{
                sw=0;
            }
            if(sw){
                for(i=0; i < ND; i++){
                    if(REAL(DATA)[i] < REAL(WHICH)[0]){
                        INTEGER(Ldata)[i]=1;
                        Lobs++;
                    }
                    else{
                        Robs++;
                    }
                }
                for(i=0; i < NTD; i++){
                    if(REAL(TDATA)[i] < REAL(WHICH)[0]){
                        INTEGER(Ltdata)[i]=1;
                    }
                }
                sob = Lobs > minbuck && Robs > minbuck;
            }
            else{
                sob=0;
            }
        }
        else{
            if(dev != 0.0 ){
                sw = (dev != 0.0) && (dev/REAL(GDEV)[0] > mindev);
                }
                else{
                    sw=0;
                }
            if(sw){
                double Lgroup[NWh];
                for(i=0; i < NWh; i++){
                    Lgroup[i]=0.0;
                }
                k=0;
                for(i=0; i < NWh; i++){
                    if(REAL(WHICH)[i] == 1.0){
                        Lgroup[k]=i+1.0;
                        k++;
                    }
                }
                for(i=0; i < k; i++){
                    for(j=0; j < ND; j++){
                        if(REAL(DATA)[j] == Lgroup[i]){
                            INTEGER(Ldata)[j]=1;
                            Lobs++;
                        }
                        else{
                            Robs++;
                        }
                    }   
                }
                for(i=0; i < k; i++){
                    for(j=0; j < NTD; j++){
                        if(REAL(TDATA)[j] == Lgroup[i]){
                            INTEGER(Ltdata)[j]=1;
                        }
                    }
                }
                sob = Lobs > minbuck && Robs > minbuck;
            }
            else{
                sob=0;
            }
        }
    }
    PROTECT(ans = allocVector(INTSXP,1));
    INTEGER(ans)[0]=sob;
    PROTECT(result_out = allocVector(VECSXP,3));
    SET_VECTOR_ELT(result_out, 0, ans);
    SET_VECTOR_ELT(result_out, 1, Ldata);
    SET_VECTOR_ELT(result_out, 2, Ltdata);
    UNPROTECT(4);
    return(result_out);  
}
