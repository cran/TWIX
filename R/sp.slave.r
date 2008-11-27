sp.slave <- function(rsp,m,test.data,Dmin=0.01,minsplit=20,
        minbucket=round(minsplit/3),topN=1,method="deviance",
        topn.method="complete",level=3,lev=0,st=1,tol=0.1,
		K=0,splitf="deviance",robust=FALSE)
{
    sp.slave.in <- function(rsp,m,test.d,dmin=Dmin,minSplit=minsplit,
                minBucket=minbucket,topn=topN,meth=method,topn.meth=topn.method,
                levelN=level,levv=lev,lstep=st,tl=tol,kfold=K,splf=splitf,rob=robust)
    {
    n <- length(m)
    E <- vector("list",3)
	Sval <-list()
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
    if(levv == 2 && meth == "grid") {
        meth <- "deviance"
        topn.meth <- "complete"
        }
    if(levv > 2 && meth == "local"){
        meth <- "deviance"
        }
    k.topn <- round(length(topn)/(n-1))
    if(k.topn == 0) k.topn <-1
	if(splf=="deviance"){
		if(rob && levv < 3){
			imp<-impor(m[,2:n],rsp,runs=10)
			TTopn<-20
			S <- lapply(m[2:n],
					splitt,
					rsp,
					meth=meth,
					lstep=lstep,
					topn=TTopn,
					topn.meth="complete",
					test=FALSE,
					level=levv,K=kfold)
			for(v in 1:length(S)){
				a<-((S[[v]]$dev/max(S[[v]]$dev)))*imp[[2]][v]*4
				S[[v]]$dev<-a
			}
		}
		else{
			if(length(topn) == 1)
				Ttest <- TRUE
			else
				Ttest <- FALSE
			S <- lapply(m[2:n],
					splitt,
					rsp,
					meth=meth,
					lstep=lstep,
					topn=if(topn.meth =="single")
							k.topn
						else
							length(topn)
					,topn.meth=topn.meth,
					test=Ttest,
					level=levv,K=kfold)
		}
	}
	if(splf == "p-adj"){
		S <- lapply(m[2:n],splitt_padj,rsp=as.numeric(rsp))
	}
	if(splf != "p-adj" & splf != "deviance"){
		stop("\n   splitf must be deviance or p-adj! \n")
	}
	if( kfold != 0 && levv < 3 && nrow(m) > 60) {
		if (S[[1]]$globD != 0 && length(S) > 0){
			k <- 0
			if (topn.meth !="single")  k <- 1
			S <- .Call("split_sum_cr",S,k,tl,PACKAGE="TWIX")
			id_Var <- S[[4]]
			k <- length(S[[1]])
			if (k != 0){
				if(k < length(topn))
					topn<-1:k
			}
		}
		else{
			S <- list()
		}
	}
	else{
		if (S[[1]]$globD != 0){
			k <- 0
			if (topn.meth !="single")  k <- 1
			S <- .Call("split_sum",S,k,tl,PACKAGE="TWIX")
			id_Var <- S[[4]]
			k <- length(S[[1]])
			if (k != 0){
				if(k < length(topn))
					topn<-1:k
			}
		}
		else{
			S <- list()
		}
	}
    splvar <- vector()
    Lcut <- Ltest <- vector("list",length(topn))
	Rcut <- Rtest <- vector("list",length(topn))
	
	make.node <- function(y,dev=S[[1]][y],gdev=S[[2]],
            spoint=S[[3]][[id_Var[y]]],var.id=S[[5]][y],
			id.var=id_Var[y],mindev=dmin,minbucket=minBucket,
			data=m,newdata=test.d,k=h,LL=levv) {
        ans <- 0
        if(length(dev) > 0 && dev > 0 && gdev > 0){
            if(is.null(dim(spoint))) {
				Splitp <- spoint[var.id]
				spvar <- data[[id.var+1]]
				Splitvar <- attr(data,"names")[id.var+1]
				ans <- .Call("split_rule",
					Splitp,
					as.character(Splitvar),
                    as.numeric(dev),
                    as.numeric(gdev),
                    as.numeric(Splitp),
                    as.integer(minbucket),
                    as.numeric(mindev),
                    as.numeric(spvar),
                    if(!is.null(newdata[[id.var+1]]))
                        as.numeric(newdata[[id.var+1]])
                    else
                        as.numeric(0.0),
					as.list(data),
					as.list(newdata),
					as.integer(LL),
                    PACKAGE="TWIX")
            }
            else{
				Splitp <- spoint[var.id,]
				spvar <- factor(data[[id.var+1]])
				Splitvar <- attr(data,"names")[id.var+1]
				attr(Splitp,"names") <- attr(spvar,"levels")
				ans <- .Call("split_rule",
					Splitp,
					as.character(Splitvar),
                    as.numeric(dev),
                    as.numeric(gdev),
                    as.numeric(Splitp),
                    as.integer(minbucket),
                    as.numeric(mindev),
                    as.numeric(spvar),
                    if(!is.null(newdata)){
                        as.numeric(factor(newdata[[id.var+1]]))
                    }else{
						as.numeric(0.0)}
					,as.list(data),
					as.list(newdata),
					as.integer(LL),
                    PACKAGE="TWIX")
            }
        } 
		else 
			ans[[1]] <- 0
        if(ans[[1]]){
            splvar[k] <<- id.var+1
            if(!is.null(newdata)) {
                Ltest[[k]]<<-ans[[5]][[1]]
                Rtest[[k]]<<-ans[[5]][[2]]
            }
            Sval[[k]] <<- ans[[6]]
            Lcut[[k]] <<- ans[[4]][[1]]
            Rcut[[k]] <<- ans[[4]][[2]]
            h <<- k+1
            return(TRUE)
        }
        else
            return(FALSE)
    }
    h <- 1
    if(length(S) > 0 && length(S[[1]]) > 0){
		for (w in topn)  {
			ans <- make.node(w)
			if(length(topn) == 1 && topn == 1 && !ans && !rob){
				w<-w+1
				if(!is.na(S[[1]][w])){
					ans<-make.node(w)
					if (!ans && !is.na(S[[1]][w+1])){
						w<-w+1
						ans<-make.node(w)
						if (!ans && !is.na(S[[1]][w+1])){
							w<-w+1
							ans<-make.node(w)
							if (!ans && !is.na(S[[1]][w+1])){
								w<-w+1
								ans<-make.node(w)
							}
						}
					}
				}
			}
		}
	}
    S <- id_Var <-0
	if(h > 1){
		E <- list(split=Sval,
            left=lapply(1:(h-1), function(z){
				if(length(Lcut[[z]][[1]]) >= minSplit && levelN > levv){
					sp.slave.in(Lcut[[z]][[1]],
                        Lcut[[z]],
                        if(length(Ltest) >= z && !is.null(Ltest[[z]])) Ltest[[z]],
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
                        kfold=kfold)
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
			right=lapply(1:(h-1),function(k) {
				if(length(Rcut[[k]][[1]]) >= minSplit && levelN > levv ){
					sp.slave.in(Rcut[[k]][[1]],
						Rcut[[k]],
						if(length(Rtest) >= k && !is.null(Rtest[[k]])) Rtest[[k]],
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
						kfold=kfold)
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
			rm(list="test.d")
		}
		else{
			Dev.test <- fit.test <- 0
		}
		rm(list="m")
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