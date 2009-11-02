predict.bundlTWIX <- function(object,newdata=NULL,sq=1:length(object$trees),ccr=FALSE,type="class", ...) {
	if(is.null(newdata))
        stop("\n newdata not supplied!\n")
    if (any(is.na(newdata)))
        stop(" missing values in newdata! \n")
	Dt <- model.frame(delete.response(terms(object$formula)),
            na.action=na.omit,newdata)
		attr(Dt,"terms") <- NULL
	Dt <- twix.data.pred(Dt,class(object)[2])
	type <- charmatch(type,c("class", "prob"))
	if (type == 1)
		pred.meth <- 0
	else if (type == 2)
		pred.meth <- 1

    k<-0
    Erg <- vector("list", length(object$trees))
	add.pred.fun <- object$add.models$predict
	pred_levels <- object[[3]][[1]][[1]][[1]]$Prob
	for(i in sq){
		add.data <- add.pred.fun(object$add.models$models[[i]],Dt)
		Bdata <- cbind(Dt,add.pred=add.data)
		ldata <- .Call("mysplit",Bdata,PACKAGE="TWIX")
		cat.levels <- attr(Dt,"cat.levels")
		mtree <- list(get.tree(object,i))
		Erg[[k<-k+1]] <- .Call("pred_TWIX",
				ldata,
				mtree,
				names(Bdata),
				cat.levels,
				length(ldata),
				as.integer(pred.meth),
				as.integer(length(pred_levels)),
				PACKAGE="TWIX")
	}
	if(pred.meth == 0){
		Erg <- as.data.frame(Erg,optional=FALSE)
		Erg <- apply(Erg,1,function(x) { names(sort(table(x),decreasing=TRUE)[1])})
	}
	else{
		Erg <- lapply(Erg,function(x,y){rownames(x[[1]])<-y;x[[1]]},y=names(pred_levels))
		name <- rownames(Erg[[1]])
		Erg <- matrix(apply(sapply(Erg,as.vector),1,sum)/length(object$trees),nrow=dim(Erg[[1]])[2],byrow=TRUE)
		colnames(Erg) <- name
	}
    if( ccr && type != 2){
        y <- newdata[names(newdata) == object$formula[[2]]]
		rate <- sum(y == Erg)/length(Erg)
        Erg <- list(predicted=Erg,CCR=rate)
        Erg
    } 
	else {
        Erg
    }
}
