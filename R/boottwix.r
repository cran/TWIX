bootTWIX<- function(formula, data=NULL, nbagg=1, topN=1, subset=NULL,
                    method="deviance", topn.method="complete", replace = TRUE, ns = 1,
                    cluster=NULL, minsplit=2, minbucket=round(minsplit/3),
                    splitf="deviance", Devmin=0.05, level=30, tol=0.01)
{
    call <- match.call()
    m <- match.call(expand=FALSE)
    m$method <- m$topn.method <- NULL
    m$cluster <- m$minsplit <- m$minbucket <- NULL
    m$Devmin <- m$topN <- m$level <- m$tol <-NULL
    m$nbagg <- m$splitf <- m$replace <- m$ns <- NULL
	
    m[[1]] <- as.name("model.frame.default")
    m <- eval(m, parent.frame())
	if(splitf == "p-adj")
        Devmin <- 1-Devmin
	test.data <- NULL
	BG <- vector("list",length=nbagg)
	if(is.null(cluster)){
		for(j in 1:nbagg){
			if(nbagg > 1)
				bag.ind <- sample(1:nrow(m), nrow(m)*ns, replace)
			else
				bag.ind <- 1:nrow(m)
			BGdata <- m[bag.ind,]
			BG[[j]] <- sp.bagg(BGdata[,1], BGdata, test.data, Devmin,
						minsplit, minbucket, cluster, topN,
						topn.method, level, method, tol, splitf)
		}
	}
	else{
		cluster_bagg <- function(n, data, t.data, devmin, Minsplit, Minbucket,
							topn, Topn.method, Level, Method, Tol, 
							Splitf, Ns, repl)
		{
			bag.ind <- sample(1:nrow(data), nrow(data)*Ns, repl)
			BGdata <- data[bag.ind,]
			BG <- sp.bagg(BGdata[,1], BGdata, test.data=NULL, Devmin=devmin,
						minsplit=Minsplit, minbucket=Minbucket, cluster=NULL, 
						topN=topn, topn.method=Topn.method, level=Level, 
						method=Method, tol=Tol, splitf=Splitf)
			BG
		}
		clusterSetupRNG(cluster)
        clusterEvalQ(cluster, library(TWIX))
		BG <- clusterApplyLB(cluster,1:nbagg,
				cluster_bagg, data=m, t.data=NULL, devmin=Devmin,
				Minsplit=minsplit, Minbucket=minbucket, topn=topN,
				Topn.method=topn.method, Level=level, Method=method,
				Tol=tol, Splitf=splitf, Ns=ns, repl=replace)
	}
	KK <- lapply(BG,function(x) x$trees)
	for(i in 1:nbagg){
		KK[[i]] <- KK[[i]][[1]]
    }
	multitree <- lapply(BG,function(y) y$multitree)
	class(KK) <- "id.tree"
	database <-list(formula=formula(terms(formula,data=m)),call=call,
			multitree=multitree,trees=KK)
	class(database) <- c("bootTWIX",splitf)
	database
}
 
