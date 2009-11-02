predict.bootTWIX <- function(object,newdata=NULL,sq=1:length(object$trees),ccr=FALSE,type="class",...){
	if(is.null(newdata))
        stop("\n newdata not supplied!\n")
	if (any(is.na(newdata)))
        stop("\n missing values in data!\n")
	Erg <- predict_boot(object,newdata,sq,ccr,type)
	if(type == "class")
		Erg <- apply(Erg,1,function(x) { names(sort(table(x),decreasing=TRUE)[1])})
	if(type == "prob"){
		name <- rownames(Erg[[1]])
		Erg <- matrix(apply(sapply(Erg,as.vector),1,sum)/length(sq),nrow=dim(Erg[[1]])[2],byrow=TRUE)
		colnames(Erg) <- name
	}
    if( ccr && type != "prob"){
        y <- newdata[names(newdata) == object$formula[[2]]]
		rate <- sum(y == Erg)/length(Erg)
        Erg <- list(predicted=Erg,CCR=rate)
        Erg
    } 
	else {
        Erg
    }
}
 



predict_boot <- function(object,newdata,sq=1:length(object$trees),ccr=FALSE,type="class", ...) {
    if (!inherits(object,"bootTWIX"))
        stop(" object not of class bootTWIX ")
	newdata <- twix.data.pred(newdata,class(object)[2])
    if (class(object)[1] == "bootTWIX") {
        Dt <- model.frame(delete.response(terms(object$formula)),
            na.action=na.omit,newdata)
    }
    else {
        Dt <- model.frame(delete.response(terms(object[[2]])),
            na.action=na.omit,newdata)
    }
	type <- charmatch(type,c("class", "prob"))
	if (type == 1)
		pred.meth <- 0
	else if (type == 2)
		pred.meth <- 1
    ##############
	fact <- sapply(Dt, function(x) !is.null(levels(x)))
	if(any(fact)){
		cat.levels <- lapply(Dt[fact], levels)
		attr(Dt, "cat.levels") <- cat.levels
	}
    ##############
    Erg <- vector()
    if (class(object)[1] == "bootTWIX") { 
        ldata <- .Call("mysplit",Dt,PACKAGE="TWIX")
		cat.levels <- attr(Dt,"cat.levels")
		mtree <- lapply(sq,function(x,y) {get.tree(y,x)},y=object)
		pred_levels <- object[[3]][[1]][[1]][[1]]$Prob
		Erg <- .Call("pred_TWIX",
				ldata,
				mtree,
				names(Dt),
				cat.levels,
				length(ldata),
				as.integer(pred.meth),
				as.integer(length(pred_levels)),
				PACKAGE="TWIX")
    }
	if(pred.meth == 0){
		colnames(Erg) <- colnames(Erg, do.NULL = FALSE, prefix = "Tree.")
		Erg <- data.frame(Erg)
	}
	else{
		Erg <- lapply(Erg,function(x,y){rownames(x)<-y;x},y=names(pred_levels))
	}
    Erg
}
