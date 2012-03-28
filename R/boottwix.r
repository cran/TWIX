bootTWIX <- function(formula, data=NULL, nbagg=25, topN=1, subset=NULL, comb=NULL,
                    method="deviance", topn.method="complete", replace = TRUE, ns = 1,
					minsplit=2, minbucket=round(minsplit/3), splitf="deviance", 
					Devmin=0.01, level=30, tol=0.01, cluster=NULL, seed.cluster=NULL)
{
    call <- match.call()
    m <- match.call(expand.dots = FALSE)
    m$method <- m$topn.method <- m$comb <- NULL
	m$minsplit <- m$minbucket <- NULL
    m$Devmin <- m$topN <- m$level <- m$tol <-NULL
    m$nbagg <- m$splitf <- m$replace <- m$ns <- NULL
	m$cluster <- m$seed.cluster <- NULL
	
	
    m[[1]] <- as.name("model.frame.default")
    m <- eval(m, parent.frame())
	if(splitf == "p-adj")
        Devmin <- 1-Devmin
	test.data <- NULL
	m <- twix.data(m,splitf)
	BG <- vector("list",length=nbagg)
	if(is.null(cluster)){
		if(is.null(comb)){
			for(j in 1:nbagg){
				if(nbagg >= 1)
					bag.ind <- sample(1:nrow(m), nrow(m)*ns, replace)
				else
					bag.ind <- 1:nrow(m)
				BGdata <- m[bag.ind,]
				BG[[j]] <- sp.bagg(BGdata[[1]], BGdata, test.data, Devmin,
						minsplit, minbucket, topN, topn.method, level, 
						method, tol, splitf)
			}
		}
		else{
			add.BG <- vector("list",length=nbagg)
			for(j in 1:nbagg){
				if(nbagg > 1)
					bag.ind <- sample(1:nrow(m), nrow(m)*ns, replace)
				else
					bag.ind <- 1:nrow(m)
				in.data <- m[bag.ind,]
				oob.data <- m[-bag.ind,]
				add.model <- comb$model(formula, data=oob.data)
				BGdata <- cbind(in.data,add.pred=comb$predict(add.model,newdata=in.data))
				add.BG[[j]] <- add.model

				BG[[j]] <- sp.bagg(BGdata[,1], BGdata, test.data, Devmin,
						minsplit, minbucket, topN, topn.method, level,
						method, tol, splitf)
			}
		}
	}
	else{
		cluster_bagg <- function(n, data, t.data, devmin, Minsplit, Minbucket,
							topn, Topn.method, Level, Method, Tol, 
							Splitf, Ns, repl, combi){
							
			bag.ind <- sample(1:nrow(data), nrow(data)*Ns, repl)
			if(is.null(combi)){
				BGdata <- data[bag.ind,]
				BG <- sp.bagg(BGdata[,1], BGdata, test.data=NULL, Devmin=devmin,
						minsplit=Minsplit, minbucket=Minbucket, 
						topN=topn, topn.method=Topn.method, level=Level, 
						method=Method, tol=Tol, splitf=Splitf)
				BG
			}
			else{
				in.data <- data[bag.ind,]
				oob.data <- data[-bag.ind,]
				add.model <- combi$model(formula, data=oob.data)
				BGdata <- cbind(in.data,add.pred=combi$predict(add.model,newdata=in.data))
				add.BG <- add.model

				BG <- sp.bagg(BGdata[,1], BGdata, test.data=NULL, Devmin=devmin,
						minsplit=Minsplit, minbucket=Minbucket,
						topN=topn, topn.method=Topn.method, level=Level, 
						method=Method, tol=Tol, splitf=Splitf)
				list(B=BG,A=add.BG)
			}
		}
		if(is.null(seed.cluster)){
			clusterSetRNGStream(cluster, iseed=sample(1:9999,6))
		}
		else{
			clusterSetRNGStream(cluster, iseed=seed.cluster)
		}
        clusterEvalQ(cluster, library(TWIX))
		BG <- clusterApplyLB(cluster,1:nbagg,
				cluster_bagg, data=m, t.data=NULL, devmin=Devmin,
				Minsplit=minsplit, Minbucket=minbucket, topn=topN,
				Topn.method=topn.method, Level=level, Method=method,
				Tol=tol, Splitf=splitf, Ns=ns, repl=replace, combi=comb)
		if(!is.null(comb)){
			add.BG <- lapply(BG,function(x) x$A)
			BG <- lapply(BG,function(x) x$B)
		}
	}
	KK <- lapply(BG,function(x) x$trees[[1]])
	multitree <- lapply(BG,function(y) y$multitree)
	class(KK) <- "id.tree"
	if(is.null(comb)){
		database <- list(formula=formula(terms(formula,data=m)),call=call,
			multitree=multitree,trees=KK)
		class(database) <- c("bootTWIX",splitf)
		database
	}
	else{
		database <- list(formula=formula(terms(formula,data=m)),call=call,
			multitree=multitree,trees=KK,add.models=list(models=add.BG,predict=comb$predict))
		class(database) <- c("bundlTWIX",splitf)
		database
	
	}
}
 








