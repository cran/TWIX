tune.TWIX <- function(formula, data = NULL, minbuck = seq(5,20,by=5), 
					maxdepth=30, Devmin = 0.01, splitf="deviance", method="deviance",
					topn.method="complete", topN=1, tol=0.15, cp=0.0, 
					xval=10, runs = 10, trace.plot = TRUE, 
					score=1, predict="best", cluster=NULL, seed.cluster=1, multicore=FALSE)
{
	call <- match.call()
	m <- match.call(expand=FALSE)
	m$xval <- m$runs <- m$data <- m$minsplit <- NULL
	m <- model.frame(formula,data)
	rsp <- model.extract(m, "response")
	n <- nrow(data)
	minbuck <- rev(minbuck)
	runs_in <- length(minbuck)
	CV <- array(0,c(runs_in,xval,runs))
	xgr <- 1:xval
	###########################################
	##### function for parallel computation of k-fold
	#####
	#####
	clus.tune.xval <- function(z, id, formula, m, minbuck, splitf, tol, cp,
							   method, topn.method, Devmin, topN, score, pr.meth)
	{
		test <- id == z
		train <- !test
		Model <- TWIX(formula, data = m[train,],level=maxdepth,minbucket=minbuck,
					minsplit=minbuck, splitf=splitf, tol=tol, method=method, score=score,
					topn.method=topn.method, Devmin=Devmin, topN=topN, cp=cp, trace=FALSE)
		if(pr.meth=="best")
			1 - predict(Model,m[test,],sq=1,ccr=TRUE)$CCR
		else
			1 - sum(as.character(bagg(Model,m[test,])) == m[test,1])/nrow(m[test,])
	}
	###########################################
	####
	for(k in 1:runs){
		id <- sample(rep(xgr, length = n), n)
		for(i in 1:runs_in) {
			if(!is.null(cluster)){
				clusterSetupRNGstream(cluster, seed=rep(seed.cluster,6))
				clusterEvalQ(cluster, library(TWIX))
				conf <- clusterApplyLB(cluster,xgr,clus.tune.xval,id,formula,m,minbuck[i],
									   splitf,tol,cp,method,topn.method,Devmin,topN,score,predict)
				conf <- unlist(conf)
			}
			else{
				if(multicore){
					conf <- mclapply(as.list(xgr),clus.tune.xval,id,formula,m,minbuck[i],
									 splitf,tol,cp,method,topn.method,Devmin,topN,score,predict,FALSE)
					conf <- unlist(conf)
				}
				else{
					conf <- rep(0,length(xgr))
					for(j in xgr) {
						test <- id == j
						train <- !test
						Model <- TWIX(formula, data = m[train,],level=maxdepth,minbucket=minbuck[i],
							minsplit=minbuck[i], splitf=splitf, tol=tol, cp=cp, method=method,
							topn.method=topn.method, Devmin=Devmin, topN=topN, score=score, trace=FALSE)
						if(predict == "best")
							conf[j] <- 1 - predict(Model,m[test,],sq=1,ccr=TRUE)$CCR
						else
							conf[j] <- 1 - sum(as.character(bagg(Model,m[test,])) == m[test,1])/sum(test)
					}
				}
			}
			CV[i,,k] <- conf
		}
	}
	CVM <- apply(CV,1,mean)
	CVVar <- apply(apply(CV,1,function(x)x),2,sd)
	if(trace.plot){
		ll <- lowess(CVM ~ minbuck)
		plot(minbuck,CVM,ylab="error",xlab="minbuck", pch=8,col=4)
		lines(ll,col=2)
	}
	out <- data.frame(minbuck=minbuck,error=CVM, sd=CVVar)
	list(best=minbuck[which.min(CVM)],tune.data=out)
}




tune.rpart <- function(formula, data = NULL, minbuck = seq(5,20,by=5), 
		parms=list(split="information"), maxdepth=30, cp = 0.0, xval=10, runs = 10,
		trace.plot = FALSE, cluster=NULL, seed.cluster=1, multicore=FALSE)
{
	call <- match.call()
	m <- match.call(expand=FALSE)
	m$xval <- m$runs <- m$data <- m$minsplit <- NULL
	m <- model.frame(formula,data)
	n <- nrow(data)
	minbuck <- rev(minbuck)
	runs_in <- length(minbuck)
	CV <- array(0,c(runs_in,xval,runs))
	xgr <- 1:xval
	###########################################
	##### parallel functions
	###########################################
	### parameter of function clus.tune.minbuck 
	### z - minbuck[i]
	### id - index-vector vor validation
	### xgr - index-vector of validation groups
	### formula - formula
	### m - data-set
	### parms - rpart parameter
	### maxdepth, cp - rpart control parameter
	###
	clus.tune.minbuck <- function(z, id, xgr, formula, m, parms, maxdepth, cp){
		conf <- rep(0,length(xgr))
		for(j in xgr){
			test <- id == j
			train <- !test
			rpt <- rpart(formula, data = m[train,],parms=parms,
				control=rpart.control(maxdepth=maxdepth, minsplit=z, minbuck=z, cp=cp, xval=0))
			conf[j] <- 1 - sum(as.character(predict(rpt,m[test,],type="class")) == as.character(m[test,1]))/sum(test)
		}
		conf
	}
	###########################################
	### parameter of function clus.tune.minbuck 
	### z - run[i]
	### id - index-vector vor validation
	### xval - number of validation groups
	### formula - formula
	### m - data-set
	### parms - rpart parameter
	### maxdepth, cp - rpart control parameter
	###
	clus.tune.runs <- function(z, xval, formula, m, minbuck, parms, maxdepth, cp){
		runs_in <- length(minbuck)
		CV <- array(0,c(runs_in,xval))
		xgr <- 1:xval
		n <- nrow(m)
		id <- sample(rep(xgr, length = n), n)
		for(i in 1:runs_in) {
			conf <- rep(0,length(xgr))
			for(j in xgr){
				test <- id == j
				train <- !test
				rpt <- rpart(formula, data = m[train,],parms=parms,
					control=rpart.control(maxdepth=maxdepth, minsplit=minbuck[i], minbuck=minbuck[i], cp=cp, xval=0))
				CV[i,j] <- 1 - sum(as.character(predict(rpt,m[test,],type="class")) == as.character(m[test,1]))/sum(test)
			}
		}
		CV
	}
	########################
	if(is.null(cluster) && !multicore){
		for(k in 1:runs){
			id <- sample(rep(xgr, length = n), n)
			for(i in 1:runs_in) {
				conf <- rep(0,length(xgr))
				for(j in xgr) {
					test <- id == j
					train <- !test
					rpt <- rpart(formula, data = m[train,],parms=parms,
						control=rpart.control(maxdepth=maxdepth,minsplit=minbuck[i],minbuck=minbuck[i],cp=cp,xval=0))
					conf[j] <- 1 - sum(as.character(predict(rpt,m[test,],type="class")) == as.character(m[test,1]))/sum(test)
				}
				CV[i,,k] <- conf
			}
		}
	}
	else if(!multicore){
			clusterSetupRNGstream(cluster, seed=rep(seed.cluster,6))
			clusterEvalQ(cluster, library(rpart))
			if(2*runs > runs_in){
				conf <- clusterApplyLB(cluster, 1:runs, clus.tune.runs, xval, formula, m, minbuck, parms, maxdepth, cp)
				for(k in 1:runs)
					CV[,,k] <- conf[[k]]
			}
			else{
				for(i in 1:runs){
					id <- sample(rep(xgr, length = n), n)
					conf <- clusterApplyLB(cluster, minbuck, clus.tune.minbuck, id, xgr, formula, m, parms, maxdepth, cp)
					for(p in 1:runs_in)
						CV[p,,i] <- conf[[p]]
				}
			}
	}
	else{
		if(2*runs > runs_in){
			conf <- mclapply(1:runs, clus.tune.runs, xval, formula, m, minbuck, parms, maxdepth, cp)
			for(k in 1:runs)
			CV[,,k] <- conf[[k]]
		}
		else{
			for(i in 1:runs){
				id <- sample(rep(xgr, length = n), n)
				conf <- mclapply(minbuck, clus.tune.minbuck, id, xgr, formula, m, parms, maxdepth, cp)
				for(p in 1:runs_in)
				CV[p,,i] <- conf[[p]]
			}
		}
	}
	CVM <- apply(CV,1,mean)
	CVVar <- apply(apply(CV,1,function(x)x),2,sd)
	if(trace.plot){
		ll <- lowess(CVM ~ minbuck)
		plot(minbuck,CVM,ylab="error",xlab="minbuck", pch=8,col=4)
		lines(ll,col=2)
	}
	out <- data.frame(minbuck=minbuck,error=CVM, sd=CVVar)
	list(best=minbuck[which.min(CVM)],tune.data=out)
}




tune.RF <- function(formula, data = NULL, mtry = 2:(ncol(data)-1), 
					ntree=500, replace=TRUE, xval=10, runs = 10, 
					trace.plot = FALSE, cluster=NULL, seed.cluster=1, 
					multicore=FALSE)
{
	call <- match.call()
	m <- match.call(expand=FALSE)
	m$xval <- m$runs <- m$data <- m$mtry <- m$replace <- m$ntree <- NULL
	m <- model.frame(formula,data)
	n <- nrow(data)
	runs_in <- length(mtry)
	CV <- array(0,c(runs_in,xval,runs))
	xgr <- 1:xval
	###########################################
	##### parallel functions
	###########################################
	### parameter of function clus.tune.mtry 
	### z - mtry[i]
	### id - index-vector vor validation
	### xgr - index-vector of validation groups
	### formula - formula
	### m - data-set
	### parms - rpart parameter
	### maxdepth, cp - rpart control parameter
	###
	clus.tune.mtry <- function(z, id, xgr, formula, m, ntree, replace){
		conf <- rep(0,length(xgr))
		for(j in xgr){
			test <- id == j
			train <- !test
			Model <- randomForest(formula, data = m[train,],mtry=z,ntree=ntree,replace=replace)
			conf[j] <- 1 - sum(predict(Model,m[test,]) == m[test,1])/sum(test)
		}
		conf
	}
	###########################################
	### parameter of function clus.tune.runs 
	### z - run[i]
	### id - index-vector vor validation
	### xval - number of validation groups
	### formula - formula
	### m - data-set
	### parms - rpart parameter
	### maxdepth, cp - rpart control parameter
	###
	clus.tune.runs <- function(z, xval, formula, m, mtry, ntree, replace){
		runs_in <- length(mtry)
		CV <- array(0,c(runs_in,xval))
		xgr <- 1:xval
		n <- nrow(m)
		id <- sample(rep(xgr, length = n), n)
		for(i in 1:runs_in) {
			conf <- rep(0,length(xgr))
			for(j in xgr){
				test <- id == j
				train <- !test
				Model <- randomForest(formula, data = m[train,],mtry=mtry[i],ntree=ntree,replace=replace)
				CV[i,j] <- 1 - sum(predict(Model,m[test,]) == m[test,1])/sum(test)
			}
		}
		CV
	}
	########################
	if(is.null(cluster) && !multicore){
		for(k in 1:runs){
			id <- sample(rep(xgr, length = n), n)
			for(i in 1:runs_in) {
				conf <- rep(0,length(xgr))
				for(j in xgr) {
					test <- id == j
					train <- !test
					Model <- randomForest(formula, data = m[train,],mtry=mtry[i],ntree=ntree)
					conf[j] <- 1 - sum(predict(Model,m[test,]) == m[test,1])/sum(test)
				}
				CV[i,,k] <- conf
			}
		}
	}
	else if(!multicore){
		clusterSetupRNGstream(cluster, seed=rep(seed.cluster,6))
		clusterEvalQ(cluster, library(randomForest))
		if(2*runs > runs_in){
			conf <- clusterApplyLB(cluster, 1:runs, clus.tune.runs, xval, formula, m, mtry, ntree, replace)
			for(k in 1:runs)
				CV[,,k] <- conf[[k]]
		}
		else{
			for(i in 1:runs){
				id <- sample(rep(xgr, length = n), n)
				conf <- clusterApplyLB(cluster, mtry, clus.tune.mtry, id, xgr, formula, m, ntree, replace)
				for(p in 1:runs_in)
					CV[p,,i] <- conf[[p]]
			}
		}
	}
	else{
		if(2*runs > runs_in){
			conf <- mclapply(1:runs, clus.tune.runs, xval, formula, m, mtry, ntree, replace)
			for(k in 1:runs)
				CV[,,k] <- conf[[k]]
		}
		else{
			for(i in 1:runs){
				id <- sample(rep(xgr, length = n), n)
				conf <- mclapply(mtry, clus.tune.mtry, id, xgr, formula, m, ntree, replace)
				for(p in 1:runs_in)
					CV[p,,i] <- conf[[p]]
			}
		}
	}
	CVM <- apply(CV,1,mean)
	CVVar <- apply(apply(CV,1,function(x)x),2,sd)
	if(trace.plot){
		ll <- lowess(CVM ~ mtry)
		plot(mtry,CVM,ylab="error",xlab="mtry", pch=8,col=4)
		lines(ll,col=2)
	}
	out <- data.frame(mtry=mtry,error=CVM, sd=CVVar)
	list(best=mtry[which.min(CVM)],tune.data=out)
}

