predict.bootTWIX <- function(object,Data,sq=1,ccr=FALSE,...) {
    Erg <- vector()
    Dt<-model.frame(delete.response(terms(object$formula)),
    na.action=na.omit,Data)
    for (i in sq){
        Erg <- cbind(Erg,apply(Dt,1,pred.value,get.tree(object,n=i)[[1]]))
        }
    if ( ccr ) {
        n <- dim(Erg)[2]
        rate <- vector(,length=n)
        y <- Data[names(Data) == object$formula[[2]]]
        for (j in 1:n) {
        rate[j] <- sum(Erg[ ,j] == y)/length(y[[1]])
        }
        colnames(Erg) <- colnames(Erg, do.NULL = FALSE, prefix = "Baum.")
        as.data.frame(Erg)
        Erg<-list(predicted=Erg,CCR=rate)
        Erg
    } else {
        colnames(Erg) <- colnames(Erg, do.NULL = FALSE, prefix = "Baum.")
        as.data.frame(Erg)
    }
}
 
