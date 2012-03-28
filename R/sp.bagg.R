
sp.bagg <- function(rsp,m,test.data=NULL,Devmin=0.01,minsplit=20,
                minbucket=10,topN=1,topn.method="complete",
                level=30, method="deviance",tol=0.0, splitf="deviance", st=1, k=0)
{
	if (method == "deviance")
		method <- 0
	else if (method == "local")
		method <- 1
	else if (method == "grid")
		method <- 2
	else
		stop("\n   method must be one of deviance, local, grid !! \n")
	if(!inherits(rsp,"factor",FALSE)) {
		stop("\n   Response must be a factor!! \n")
	}
		
	sp.twix <- function(rsp,m, test.d=NULL, dmin=Devmin, minSplit=minsplit,
                minBucket=minbucket, topn=topN, topn.meth=topn.method,
                levelN=level, meth=method, tl=tol, splf=splitf, lstep=st, K=k,
				lev=0)
	{
    n <- length(m)
    E <- vector("list",3)
    lev <- lev+1
    if (length(topN) == 1) {
        topn <- 1:topN
    }
    else {
        if (!is.na(topN[lev]) && topN[lev] != 1) {
            topn <- 1:topN[lev]
        }
        else
            topn <- 1
    }
    if(lev == 2 && meth == 2){
        meth <- 0
        #topn.meth <- "complete"
        }
    if(lev > 2 && meth == 1){
        meth <- 0
        }
	n_cut <- n_topn <- length(topn)
	if(topn.meth == "single"){
		n_cut <- round(n_topn/(n-1))
		if(n_cut == 0){
			n_cut <- 1
		}
	}
	if(splf == "deviance"){
		if(meth == 2 || (K != 0 && lev < 3)){
			S <- lapply(m[2:n],
						splitt_cross,
						rsp,
						meth=meth,
						lstep=lstep,
						topn=n_cut,
						topn.meth=topn.meth,
						test=FALSE,
						level=lev,K=K)
		}
		else{
			S <- .Call("var_split_dev", m, as.integer(meth),
						as.integer(lstep), as.integer(n_cut),
						as.logical(FALSE), as.numeric(K),
						as.integer(minBucket), as.integer(lev),
						as.logical(TRUE),as.numeric(tl),
						PACKAGE="TWIX")
		}
		S <- .Call("split_summary_dev", S, as.numeric(tl), PACKAGE="TWIX")
	}
	if(splf == "p-adj"){
		S <- .Call("var_split_adj", m, as.numeric(0.1), as.numeric(0.9),
				   as.logical(FALSE), as.logical(TRUE), as.numeric(tl),
				   as.integer(minBucket), PACKAGE="TWIX")
		S <- .Call("split_summary_padj", S, as.numeric(tl), PACKAGE="TWIX")
	}
	k <- length(S[[1]])
	if(k < length(topn)){
		topn <- 1:k
	}
	Sval <- Rcut <- Rtest <- Lcut <- Ltest <- list()
    h_level <- 1
	erg_node <- .Call("make_node", topn, S, as.integer(minBucket), as.numeric(dmin),
			m, test.d, as.integer(lev), environment(), PACKAGE="TWIX")
	S <- 0
	if (h_level > 1) {
		E <- list(split=Sval,
		left=my.lapply(1:(h_level-1),
				function(z) {
				if(length(Lcut[[z]][[1]]) >= minSplit && levelN > lev){
				sp.twix(Lcut[[z]][[1]],
					Lcut[[z]],
					if(length(Ltest) >= z && length(Ltest[[z]][[1]]) > 0) Ltest[[z]],
					dmin=dmin,
					minSplit=minSplit,
					minBucket=minBucket,
					topn=topn,
					meth=meth,
					topn.meth=topn.meth,
					levelN=levelN,
					lev=lev,
					lstep=lstep,
					tl=tl,
					K=K)
				}
				else {
					yLcut <- Lcut[[z]][[1]]
					TB.L <- .Call("tw_table",yLcut,PACKAGE="TWIX")
					Pred.class <- attr(which.max(TB.L),"names")
					if(FALSE){#if(length(Ltest) > z){
						yLtest <- Ltest[[z]][[1]]
						Dev.test <- .Call("Dev_leaf",yLtest,PACKAGE="TWIX")
						id.test <- Pred.class == yLtest
						fit.test <- .Call("tw_table",yLtest[id.test],PACKAGE="TWIX")
					}
					else {
						Dev.test <- fit.test <- 0
					}
					list(Obs=sum(TB.L),
						Prob=TB.L/sum(TB.L),
						Pred.class=Pred.class,
						Dev=.Call("Dev_leaf",yLcut,PACKAGE="TWIX"),
						Dev.test=Dev.test,
						fit.tr=sum(Pred.class == yLcut),
						fit.test=fit.test)
				}
			}
		),
        right=my.lapply(1:(h_level-1),
            function(k) {
            if(length(Rcut[[k]][[1]]) > minSplit && levelN > lev){
                sp.twix(Rcut[[k]][[1]],
                Rcut[[k]],
                if(length(Rtest) > k && length(Rtest[[k]][[1]]) > 0) Rtest[[k]],
                        dmin=dmin,
                        minSplit=minSplit,
                        minBucket=minBucket,
                        topn=topn,
                        meth=meth,
                        topn.meth=topn.meth,
                        levelN=levelN,
                        lev=lev,
                        lstep=lstep,
                        tl=tl,
                        K=K)
            }
			else{
				yRcut <- Rcut[[k]][[1]]
                TB.R <- .Call("tw_table",yRcut,PACKAGE="TWIX")
                Pred.class <- attr(which.max(TB.R),"names")
                if(FALSE){#if(length(Rtest) > k){
					yRtest <- Rtest[[k]][[1]]
                    Dev.test <- .Call("Dev_leaf",yRtest,PACKAGE="TWIX")
                    id.test <- Pred.class == yRtest
                    fit.test <- .Call("tw_table",yRtest[id.test],PACKAGE="TWIX")
                }
				else{
                    Dev.test <- fit.test <- 0
                }
                list(Obs=sum(TB.R),
                    Prob=TB.R/sum(TB.R),
                    Pred.class=Pred.class,Dev=.Call("Dev_leaf",yRcut,PACKAGE="TWIX"),
                    Dev.test=Dev.test,
                    fit.tr=sum(Pred.class == yRcut),
                    fit.test=fit.test)
              }
             }
          )
        )
    } else {
		ylev <- .Call("tw_table",rsp,PACKAGE="TWIX")
		Pred.class <- attr(which.max(ylev),"names")
		Dev <- .Call("Dev_leaf",rsp,PACKAGE="TWIX")
		fit.tr <- sum(Pred.class == rsp)
		if(FALSE){#if(length(test.d[[1]]) > 0) {
			ytest <- test.d[[1]]
			Dev.test <- .Call("Dev_leaf",ytest,PACKAGE="TWIX")
			id.test <- Pred.class == ytest
			fit.test <- .Call("tw_table",ytest[id.test],PACKAGE="TWIX")
		} else {
			Dev.test <- fit.test <- 0
			}
		m <- test.d <- ytest <- 0
		E <- list(Obs=sum(ylev),Prob=ylev/sum(ylev),
					Pred.class=Pred.class,Dev=Dev,
					Dev.test=Dev.test,fit.tr=fit.tr,
					fit.test=fit.test)
		}
	E
	}

	K <- sp.twix(rsp,m)
    if(!is.null(K$split)){
		if(splitf != "p-adj"){
			KK <- s_bag(K)
		}
		else{
			KK <- s_bag_padj(K)
		}
        id.leng <- sapply(KK,function(x) length(x$id))
        id.out <- id.leng == 1
        if(sum(id.out) > 0 && sum(id.out) < length(KK)){
            KK <- KK[!id.out]
            id.leng <- id.leng[!id.out]
        }
        if(sum(id.out) == length(id.leng)){
            KK <- list(KK[[1]])
        }
		ntr <- dim(m)[1]
		ntest <- 1
        for(i in 1:length(KK)) {
            KK[[i]]$dev.test <- K$split[[1]]$Dev.test-KK[[i]]$dev.test
            KK[[i]]$fit.tr <- sum(KK[[i]]$fit.tr)/ntr
            KK[[i]]$fit.test <- sum(KK[[i]]$fit.test)/ntest
        }
    }
    else{
		l <- list(dev.test=0, dev=0, id=0)
        K <- c(l,K)
        KK <- list(K)
        K <- list(split=list(K))
    }
	if(length(KK) > 1){
		KK <- KK[my.sort(sapply(KK,function(x)x$dev),
			index.return=TRUE,decreasing=TRUE)$ix]
		KK <- KK[my.sort(sapply(KK,function(x)x$fit.tr),
			index.return=TRUE,decreasing=TRUE)$ix]
	}
	
    database <- list(multitree=K,trees=KK)
	if(splitf != "p-adj")
		class(database) <- c("TWIX","deviance")
	else
		class(database) <- c("TWIX","p-adj")
    database
}

