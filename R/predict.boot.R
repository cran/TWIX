predict_boot <- function(object,newdata,sq=1,ccr=FALSE,type="class", ...) {
    if (!inherits(object,"bootTWIX"))
        stop(" object not of class bootTWIX ")
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
    ############## make 'data.frame' numeric 
    fact <- sapply(Dt, function(x) !is.null(levels(x)))
    char <- sapply(Dt, is.character)
    sDt <- 1:ncol(Dt)
    if(any(fact | char)) {
        for (j in sDt[char])
            Dt[[j]] <- as.factor(Dt[[j]])
        fact <- fact | char
        cat.levels <- lapply(Dt[fact], levels)
        for (j in sDt[fact])
            Dt[[j]] <- as.integer(Dt[[j]])
        attr(Dt, "cat.levels") <- cat.levels
    }
    ##############
    Erg <- vector()
    if (class(object)[1] == "bootTWIX") { 
        ldata <- split(Dt,1:nrow(Dt))
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
    if(ccr){
        n <- dim(Erg)[2]
        rate <- vector(,length=n)
        y <- newdata[names(newdata) == object$formula[[2]]]
        for (j in 1:n) {
            rate[j] <- sum(Erg[ ,j] == y[[1]])/length(y[[1]])
            }
		if(pred.meth == 0){
			colnames(Erg) <- colnames(Erg, do.NULL = FALSE, prefix = "Tree.")
			Erg <- data.frame(Erg)
			Erg <- list(predicted=Erg,CCR=rate)
		}
		else{
			Erg <- lapply(Erg,function(x,y){rownames(x)<-y;x},y=names(pred_levels))
		}
    }
    else {
		if(pred.meth == 0){
			colnames(Erg) <- colnames(Erg, do.NULL = FALSE, prefix = "Tree.")
			Erg <- data.frame(Erg)
		}
		else{
			Erg <- lapply(Erg,function(x,y){rownames(x)<-y;x},y=names(pred_levels))
		}	
	}
    Erg
}
