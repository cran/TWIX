predict.bootTWIX <- function(object,newdata=NULL,sq=1,ccr=FALSE,...){
	if(is.null(newdata))
        stop("\n newdata not supplied!\n")
	if (any(is.na(newdata)))
        stop("\n missing values in data!\n")
	Erg <- predict_boot(object,newdata,sq)
	Erg <- apply(Erg,1,function(x) { names(sort(table(x),decreasing=TRUE)[1])})
    if( ccr ){
        y <- newdata[names(newdata) == object$formula[[2]]]
		rate <- sum(y == Erg)/length(Erg)
        Erg <- list(predicted=Erg,CCR=rate)
        Erg
    } 
	else {
        Erg
    }
}
 
