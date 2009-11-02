sp.slave <- function(rsp, m, test.data, Dmin=0.01, minsplit=20,
        minbucket=round(minsplit/3), topN=1, method=0,
        topn.method="complete", level=30, lev=0, st=1, tol=0.15,
		K=0, splitf="deviance", robust=FALSE)
{
    sp.slave.in <- function(rsp,m,test.d,dmin=Dmin,minSplit=minsplit,
                minBucket=minbucket,topn=topN,meth=method,topn.meth=topn.method,
                levelN=level,levv=lev,lstep=st,tl=tol,KK=K,splf=splitf,rob=robust)
    {
    n <- length(m)
    E <- vector("list",3)
    levv <- levv+1
    if(length(topN) == 1) {
        topn <- 1:topN
    }
    else {
        if (!is.na(topN[levv]) && topN[levv] != 1) {
            topn <- 1:topN[levv]
        }
        else
            topn <- 1
    }
    if(levv == 2 && meth == 2) {
        meth <- 0
        #topn.meth <- "complete"
        }
    if(levv > 2 && meth == 1){
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
		if(rob && levv < 3){
			imp <- impor(m[,2:n],rsp,runs=10)
			S <- .Call("var_split_dev", m, as.integer(meth), 
					   as.integer(lstep), as.integer(n_cut),
					   as.logical(FALSE), as.numeric(KK),
					   as.integer(minBucket), as.integer(lev),
					   PACKAGE="TWIX")
			for(v in 1:length(S)){
				a <- ((S[[v]][[1]]/max(S[[v]][[1]])))*imp[[2]][v]*4
				a[is.na(a)] <- 0
				S[[v]][[1]] <- a
			}
		}
		else{
			S <- .Call("var_split_dev", m, as.integer(meth), 
					   as.integer(lstep), as.integer(n_cut),
					   as.logical(FALSE), as.numeric(KK),
					   as.integer(minBucket), as.integer(lev),
					   PACKAGE="TWIX")
		}
	}
	else{
		S <- .Call("var_split_adj", m, as.numeric(0.1), as.numeric(0.9),
					as.logical(FALSE), PACKAGE="TWIX")
	}
	S <- .Call("split_sum",S,as.numeric(tl),PACKAGE="TWIX")
	k <- length(S[[1]])
	if(k < length(topn)){
		topn <- 1:k
	}
	Sval <- Rcut <- Rtest <- Lcut <- Ltest <- list()
	h_level <- 1
	erg_node <- .Call("make_node", topn, S, as.integer(minBucket), as.numeric(dmin),
					m, test.d, as.integer(lev), environment(), PACKAGE="TWIX")

	if(h_level > 1){
		E <- list(split=Sval,
            left=my.lapply(1:(h_level-1), function(z){
				if(length(Lcut[[z]][[1]]) >= minSplit && levelN > levv){
					sp.slave.in(Lcut[[z]][[1]],
                        Lcut[[z]],
                        if(length(Ltest) >= z && length(Ltest[[z]][[1]]) > 0) Ltest[[z]],
                        dmin=dmin,
                        minSplit=minSplit,
                        minBucket=minBucket,
                        topn=topn,
                        meth=meth,
                        topn.meth=topn.meth,
                        levelN=levelN,
                        levv=levv,
                        lstep=lstep,
                        tl=tl,
                        KK=KK)
                    }
					else{
						yLcut <- Lcut[[z]][[1]]
                        TB.L <- .Call("tw_table",yLcut,PACKAGE="TWIX")
                        Pred.class <- attr(which.max(TB.L),"names")
                        if(length(Ltest) >= z && !is.null(Ltest[[z]][[1]])){
							yLtest <- Ltest[[z]][[1]]
                            Dev.test <- .Call("Dev_leaf",yLtest,PACKAGE="TWIX")
                            id.test <- Pred.class == yLtest
                            fit.test <- sum(id.test)
                        }
						else{
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
			right=my.lapply(1:(h_level-1),function(k) {
				if(length(Rcut[[k]][[1]]) >= minSplit && levelN > levv ){
					sp.slave.in(Rcut[[k]][[1]],
						Rcut[[k]],
						if(length(Rtest) >= k && length(Rtest[[k]][[1]]) > 0) Rtest[[k]],
						dmin=dmin,
						minSplit=minSplit,
						minBucket=minBucket,
						topn=topn,
						meth=meth,
						topn.meth=topn.meth,
						levelN=levelN,
						levv=levv,
						lstep=lstep,
						tl=tl,
						KK=KK)
				} 
				else{
					yRcut <- Rcut[[k]][[1]]
					TB.R <- .Call("tw_table",yRcut,PACKAGE="TWIX")
					Pred.class <- attr(which.max(TB.R),"names")
					if(length(Rtest) >= k && !is.null(Rtest[[k]])){
						yRtest <- Rtest[[k]][[1]]
						Dev.test <- .Call("Dev_leaf",yRtest,PACKAGE="TWIX")
						id.test <- Pred.class == yRtest
						fit.test <- sum(id.test)
					}
					else{
						Dev.test <- fit.test <- 0
					}
					list(Obs=sum(TB.R),
						Prob=TB.R/sum(TB.R),
						Pred.class=Pred.class,
						Dev=.Call("Dev_leaf",yRcut,PACKAGE="TWIX"),
						Dev.test=Dev.test,
						fit.tr=sum(Pred.class == yRcut),
						fit.test=fit.test)
				}
			}
          )
        )
	}
	else{
		ylev <- .Call("tw_table",rsp,PACKAGE="TWIX")
		Pred.class <- attr(which.max(ylev),"names")
		Dev <- .Call("Dev_leaf",rsp,PACKAGE="TWIX")
		fit.tr <- sum(Pred.class == rsp)
		if(length(test.d[[1]]) > 0) {
			ytest <- test.d[[1]]
			Dev.test <- .Call("Dev_leaf",ytest,PACKAGE="TWIX")
			fit.test <- sum(Pred.class == ytest)
		}
		else{
			Dev.test <- fit.test <- 0
		}
		E <- list(Obs=sum(ylev),
				Prob=ylev/sum(ylev),
				Pred.class=Pred.class,
				Dev=Dev,
				Dev.test=Dev.test,
				fit.tr=fit.tr,
				fit.test=fit.test)
	}
	E
	}
	sp.slave.in(rsp,m,test.data)
}