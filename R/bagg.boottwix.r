bagg <- function(object, ...) UseMethod("bagg")

bagg.bootTWIX <- bagg.TWIX <-
    function(object, ...) bagg.default(object, ...)

bagg.default <- function(object, newdata=NULL, sq=1:length(object$trees), 
							aggregation="weighted", type="class", ...)
{
    if(is.null(newdata))
        stop("Data not supplied")
    if (any(is.na(newdata)))
        stop("\n missing values in data!\n")
	if (!is.data.frame(newdata)) 
		newdata <- as.data.frame(newdata)
    if(!inherits(object, "TWIX") & !inherits(object, "bootTWIX"))
        stop("Not legitimate object")
	if(type != "class" & type != "prob")
		stop("\n   'type' must be 'class' or 'prob'! \n")
	if(inherits(object, "bootTWIX")){
		tpred <- predict_boot(object,newdata,sq)
		out <- apply(tpred,1,function(x){names(sort(table(x),decreasing=TRUE)[1])})
	}
	else{
		agg.type <- charmatch(aggregation,c("majority", "weighted", "BMA"))
		out <- comb.method(agg.type, object, newdata, sq, type)
	}
	out
}






comb.method <- function(agg.type, object, newdata, sq, type)
{
### agg.type = 1  - majority voting
###
	if(agg.type == 1){
		tpred0 <- predict(object, newdata, sq, FALSE, type)
		if(type == "class"){
			out <- apply(tpred0,1,function(x){names(sort(table(x),decreasing=TRUE)[1])})
		}
		if(type == "prob"){
			if(length(sq) > 1){
				out <- matrix(apply(matrix(unlist(tpred0,use.names=FALSE),ncol=length(tpred0)),1,sum),ncol=dim(tpred0[[1]])[2])
				out <- out/length(tpred0)
				colnames(out) <- dimnames(tpred0[[1]])[[2]]
			}else{
				out <- tpred0
			}
		}
	}
### agg.type = 2  - Bayesian combination 1
###
	if(agg.type == 2){
		certainty <- sapply(object$trees,function(x) prod(x$y.prob[x$Obs!=0]))
		tpred <- predict(object, newdata, sq, type="prob")
		if(length(sq) > 1){
			out <- matrix(apply(mapply(get("*"), certainty[sq],tpred),1,sum),nrow=nrow(newdata))
			out <- out/length(tpred)
			if(type == "prob"){
				colnames(out) <- dimnames(tpred[[1]])[[2]]
			}
			if(type == "class"){
				out <- (dimnames(tpred[[1]])[[2]])[apply(out,1,which.max)]
			}
		}else{
			out <- certainty[sq]*tpred
			if(type == "class"){
				out <- (dimnames(tpred)[[2]])[apply(out,1,which.max)]
			}
		}
	}
### agg.type = 3  -  Dempster-Shafer combination 
###
	if(agg.type == 3){
		certainty <- sapply(object$trees,function(x) prod(x$y.prob[x$Obs!=0]))
		tpred <- predict(object, newdata, sq, type="prob")
		
		answ <- certainty[1]*tpred[[1]]
		for(p in sq[-1]){
			answ <- answ + certainty[p]*tpred[[p]]
		}
		out <- (dimnames(tpred[[1]])[[1]])[apply(answ,2,which.max)]
	}
	if(is.na(agg.type)){
		stop("\n aggregation must be one of 'weighted', 'majority' or 'BMA'")
	}
	out
}