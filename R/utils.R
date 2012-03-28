nlogn <- function(x) { if (x == 0) 0 else x*log(x) }

fullrks <- function(m) {
	if (!is.factor(m)) sort(m, method="quick", index.return=TRUE)[[2]] else 0
}


my.sort <- function(x,index.return=FALSE,decreasing=FALSE) {
	if(decreasing) 
		x <- -x
	y <- sort(x, method="quick", index.return=index.return)
	if(decreasing) 
		y$x <- -y$x
	y
}


my.lapply <- function (X, FUN, ...) 
{
    FUN <- match.fun(FUN)
	class(X) <- "list"
    lapply(X, FUN)
}


my.design <- function(x){
	lev <- levels(x)
	if(length(lev) > 2){
		ans <- sapply(lev,function(x,y) as.numeric(x==y),y=x)
		as.data.frame(ans)
	}
	else{
		as.numeric(x)
	}
}


twix.data <- function(Dt,splitf){
	fact <- sapply(Dt, function(x) !is.null(levels(x)))
	char <- sapply(Dt, is.character)
	int <- sapply(Dt, is.integer)
	sDt <- 1:ncol(Dt)
	if(any(fact | char)) {
		for (j in sDt[char])
			Dt[[j]] <- factor(Dt[[j]])
	}
	if(any(int)) {
		for (j in sDt[int])
			Dt[[j]] <- as.numeric(Dt[[j]])
	}
	if(splitf == "p-adj"){
		if(sum(fact) > 1){
			fact[1] <- FALSE
			Dt_dummy <- lapply(Dt[fact],my.design)
			Dt <- cbind(Dt[!fact],as.data.frame(Dt_dummy))
		}
	}
	Dt
}



twix.data.pred <- function(newdata,splitf){ 
	fact <- sapply(newdata, function(x) !is.null(levels(x)))
	char <- sapply(newdata, is.character)
	int <- sapply(newdata, is.integer)
	sDt <- 1:ncol(newdata)
	if(any(fact | char)) {
		for (j in sDt[char])
			newdata[[j]] <- as.factor(newdata[[j]])
		fact <- fact | char
		if(splitf == "p-adj"){
			Dt_dummy <- lapply(newdata[fact],my.design)
			newdata <- cbind(newdata[!fact],as.data.frame(Dt_dummy))
		}
		cat.levels <- lapply(newdata[fact], levels)
		attr(newdata, "cat.levels") <- cat.levels
	}
	if(any(int)) {
		for (j in sDt[int])
			newdata[[j]] <- as.numeric(newdata[[j]])
	}
	newdata
}











