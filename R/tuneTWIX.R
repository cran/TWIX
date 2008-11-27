tuneTWIX <- function(formula, data = NULL, minbuck = seq(5,20,by=5), 
					maxdepth=30, xval = 10, runs = 10, trace.plot = TRUE)
{
	call <- match.call()
	m <- match.call(expand=FALSE)
	m$xval <- m$runs <- m$data <- m$minsplit <- NULL
	m <- model.frame(formula,data)
	rsp <- model.extract(m, "response")
	n <- nrow(data)
	minbuck <- rev(minbuck)
	runs_in <- length(minbuck)
	CV <- array(0,c(runs,runs_in))
	xgr <- 1:xval
	####
	for(k in 1:runs){
		id <- sample(rep(xgr, length = n), n)
		for(i in 1:runs_in) {
			conf <- rep(0,length(xgr))
			for(j in xgr) {
				test <- id == j
				train <- !test
				rpt <- rpart(formula, data=m[train,],parms=list(split="information"),
					control=rpart.control(cp=0.0001,maxdepth=maxdepth,minbucket=minbuck[i],
					minsplit=minbuck[i],xval=10))
				conf[j] <- sum(predict(rpt,newdata=m[test,], type="class") != m[test,1])/sum(test)
			}
			CV[k,i] <- mean(conf)
		}
	}
	CV <- apply(CV,2,mean)
	ll <- lowess(CV ~ minbuck)
	if(trace.plot){
		plot(minbuck,CV,ylab="error",xlab="minbuck", pch=8,col=4)
		lines(ll,col=2)
	}
	out <- data.frame(minbuck=minbuck[which.min(CV)],error=min(CV))
	out
}
