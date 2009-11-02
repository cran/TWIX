bagg <- function(object, ...) UseMethod("bagg")

bagg.bootTWIX <- bagg.TWIX <-
    function(object, ...) bagg.default(object, ...)

bagg.default <- function(object,newdata=NULL,sq=1:length(object$trees),aggregation="weighted", ...) {
    if(is.null(newdata))
        stop("Data not supplied")
    if (any(is.na(newdata)))
        stop("\n missing values in data!\n")
	if (!is.data.frame(newdata)) 
		newdata <- as.data.frame(newdata)
    if(!inherits(object, "TWIX") & !inherits(object, "bootTWIX"))
        stop("Not legitimate object")
	if(inherits(object, "bootTWIX")){
		tpred <- predict_boot(object,newdata,sq)
		out <- apply(tpred,1,function(x) { names(sort(table(x),decreasing=TRUE)[1])})
	}
	else{
		agg.type <- charmatch(aggregation,c("weighted", "majority"))
		if(agg.type == 1){
			certainty <- sapply(object$trees,function(x) x$fit.tr*(prod(x$sd.tr[x$Obs!=0])))
			tpred <- predict(object,newdata,sq,type="prob")
			
			answ <- certainty[1]*tpred[[1]]
			for(p in sq[-1]){
				answ <- answ + certainty[p]*tpred[[p]]
			}
			out <- (dimnames(tpred[[1]])[[1]])[apply(answ,2,which.max)]
		}
		if(agg.type == 2){
			tpred <- predict(object,newdata,sq)
			out <- apply(tpred,1,function(x) { names(sort(table(x),decreasing=TRUE)[1])})
		}
		if(is.na(agg.type)){
			stop("\n aggregation must be one of 'weighted', 'majority'")
		}
	}
	out
}
