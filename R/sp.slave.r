sp.slave <- function(rsp,m,test.data,Dmin=0.01,minsplit=20,
        minbucket=round(minsplit/3),topN=1,method="deviance",
        topn.method="complete",level=3,lev=0,st=1,tol=0.1,K=0,oldspvar=0)
{
	n <- length(m)
    sp.sl <- function(rsp,m,test,dmin=Dmin,minSplit=minsplit,minBucket=minbucket,
            topn=topN,meth=method,topnmeth=topn.method,
            levelN=level,levv=lev,lstep=st,tl=tol,kfold=K,oldspV=oldspvar)
    {
    Dev.leaf <- function(x) {
        TD <- 0
        CCR <- .Call("tw_table",as.integer(x),levels(x),PACKAGE="TWIX")
        s <- sum(CCR)
        TD <- -sum(s*sapply(CCR,function(x,y){ if(x/y == 0) 0 else (x/y)*log(x/y)},y=s))
        round(TD,digits=6)
    }
    E <- Sval <- list()
    levv <- levv+1
    if (length(topN) == 1) {
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
        topnmeth <- "complete"
        }
    k.topn <- round(length(topn)/(n-1))
    if(k.topn == 0) k.topn <-1

	S <- lapply(m[2:n],
		splitt,
		rsp,
        meth=meth,
        lstep=lstep,
        topn=if(topnmeth =="single")
        			k.topn
            	else
                    length(topn)
		,topn.meth=topnmeth,
		test=FALSE,
		level=levv,K=kfold)


    if( kfold != 0 && levv < 2 && nrow(m) > 60) {
        S_summ <- Var_id <- id_Var <- vector()
        which <- list()
        globD <-k<-0
        if (S[[1]]$globD != 0){
            if (topnmeth !="single")  k <- 1
            S_summ <- .Call("split_sum_cr",S,length(S),as.integer(k),as.numeric(tl),PACKAGE="TWIX")
            S <- 0
            which <- S_summ[[3]]
            Var_id <- S_summ[[5]]
            id_Var <- S_summ[[4]]
            S_summ <- S_summ[1:2]
            k <- length(S_summ[[1]])
            if (k != 0){
                if(k < length(topn))
                    topn<-1:k
                }
            }
    }
    else {
        Var_id <- id_Var <- vector()
        S_summ <- which <- list()
        k<-0
        if (S[[1]]$globD != 0){
            if (topnmeth !="single")  k <- 1
            S_summ <- .Call("split_sum",S,length(S),as.integer(k),as.numeric(tl),PACKAGE="TWIX")
            S <- 0
            which <- S_summ[[3]]
            Var_id <- S_summ[[5]]
            id_Var <- S_summ[[4]]
            S_summ <- S_summ[1:2]
            k <- length(S_summ[[1]])
            if (k != 0){
                if(k < length(topn))
                        topn<-1:k
			}
		}
	}
print("HALLO0")
    splvar <- vector()
    Lcut <- Ltest <- Rcut <- Rtest <- list()
    make.node<- function(y,dev=S_summ[[1]][y],gdev=S_summ[[2]],
            spoint=which[[id_Var[y]]],var.id=Var_id[y],id.var=id_Var[y],
            mindev=dmin,minbucket=minBucket,data=m,newdata=test,k=h,oldspv=oldspV,LL=levv) {
        ans <- 0
        if(LL <= 2) gdev <- dev
        if(dim(data)[2] <= 2) oldspv <- 0
        if(id.var != oldspv){
            if(is.null(dim(spoint))) {
            ans <- .Call("split_rule",
                    as.numeric(dev),
                    as.numeric(gdev),
                    as.numeric(spoint[var.id]),
                    as.integer(minbucket),
                    as.numeric(mindev),
                    as.numeric(data[,id.var+1]),
                    as.numeric(newdata[,id.var+1]),
                    PACKAGE="TWIX")
            Splitp <- spoint[var.id]
            Splitvar <- attr(data[id.var+1],"names")
            spvar <- data[,id.var+1]
            dist <- quantile(c(Splitp-spvar[spvar < Splitp],spvar[spvar > Splitp]-Splitp))[3]
            }
            else{
            ans <- .Call("split_rule",
                    as.numeric(dev),
                    as.numeric(gdev),
                    as.numeric(spoint[var.id,]),
                    as.integer(minbucket),
                    as.numeric(mindev),
                    as.numeric(factor(data[,id.var+1])),
                    if(!is.null(newdata)){
                        as.numeric(factor(newdata[,id.var+1]))
                    }else{
                    as.numeric(0.0)}
                    ,PACKAGE="TWIX")
            Splitp <- spoint[var.id,]
            Splitvar <- attr(data[id.var+1],"names")
            attr(Splitp,"names")<-levels(factor(data[[id.var+1]]))
            dist<-ks.t<-0
            }
        } else ans[[1]] <- 0
        if(ans[[1]]){
            rsp <- data[[1]]
            ylev <- .Call("tw_table",as.integer(rsp),levels(rsp),PACKAGE="TWIX")
            Prob <- ylev/nrow(data)
            Pred.class <-attr(which.max(ylev),"names")
            if(is.null(dim(spoint))){
                max.cl.l<-names(sort(table(rsp[spvar < Splitp]),decreasing=TRUE)[1])
                max.cl.r<-names(sort(table(rsp[spvar >= Splitp]),decreasing=TRUE)[1])
                id.p.l <- which(spvar < Splitp & rsp == max.cl.l)
                id.p.r <- which(spvar >= Splitp & rsp != max.cl.l & rsp == max.cl.r)
                if(length(id.p.l) > 10 && length(id.p.r) > 10){
                    ks.t <- mean(spvar[id.p.r])-mean(spvar[id.p.l])
                }else{
                    ks.t<-0
                }
            }
            splvar[k] <<- id.var+1
            if(nrow(newdata) > 0 && !is.null(newdata)) {
                rspt <- newdata[[1]]
                dev.test <- Dev.leaf(rspt)
                id.test <- Pred.class == rspt
                fit.test <- .Call("tw_table",as.integer(rspt[id.test]),levels(rspt[id.test]),PACKAGE="TWIX")
                rspt <- 0
                Ltest[[k]]<<-newdata[as.logical(ans[[3]]),]
                Rtest[[k]]<<-newdata[!as.logical(ans[[3]]),]
            }
            else {
                dev.test <- fit.test <- 0
                Ltest[[k]] <<- NULL
                Rtest[[k]] <<- NULL
                }
            Sval[[k]] <<- list(Splitp=Splitp,Splitvar=Splitvar,
                    Obs=sum(ylev),Dev=dev,
                    Dev.test= dev.test,
                    fit.tr=sum(rep(Pred.class,length(rsp)) == rsp),
                    fit.test=fit.test,TD=S_summ[[2]],
                    Pred.class=Pred.class,Prob=Prob,dist=dist,ks.t=ks.t)
            Lcut[[k]] <<- data[as.logical(ans[[2]]),]
            Rcut[[k]] <<- data[!as.logical(ans[[2]]),]
            h <<- k+1
            return(TRUE)
        }
        else
            return(FALSE)
    }
print("HALLO1")
    h<-1
    if(length(S_summ) > 0 && length(S_summ[[1]]) > 0)
    for (w in topn)  {
        ans <- make.node(w)
        if(length(topn) == 1 && topn == 1 && !ans){
            w<-w+1
            if(!is.na(S_summ[[1]][w])){
                ans<-make.node(w)
                if (!ans && !is.na(S_summ[[1]][w+1])){
                    w<-w+1
                    make.node(w)
                }
            }
        }
    }
print("HALLO2")
    S_summ<-which<-Var_id<-id_Var<-0
	if (length(Lcut) > 0 && length(Rcut) > 0 ) {
		m<-test<-0
		E <-list(split=Sval,
            left=lapply(1:length(Lcut),function(z) {
				if (length(Lcut[[z]][[1]]) > minSplit && levelN > levv)
				{
					sp.sl(Lcut[[z]][[1]],Lcut[[z]],
                        if(length(Ltest) > 0) Ltest[[z]],
                        dmin=dmin,
                        minSplit=minSplit,
                        minBucket=minBucket,
                        topn=topn,
                        meth=meth,
                        topnmeth=topnmeth,
                        levelN=levelN,
                        levv=levv,
                        lstep=lstep,
                        tl=tl,
                        kfold=kfold,oldspV=splvar[z])
				}
				else {
					TB.L <- .Call("tw_table",as.integer(Lcut[[z]][[1]]),levels(Lcut[[z]][[1]]),PACKAGE="TWIX")
                        Pred.class <-attr(which.max(TB.L),"names")
                        if(length(Ltest) > 0 && length(Ltest[[z]][[1]]) > 0){
                            Dev.test <-Dev.leaf(Ltest[[z]][[1]])
                            id.test <- Pred.class == Ltest[[z]][[1]]
                            fit.test <- .Call("tw_table",as.integer(Ltest[[z]][[1]][id.test]),
                                levels(Ltest[[z]][[1]][id.test]),PACKAGE="TWIX")
                        } else {
                            Dev.test <- fit.test <-0
                        }
                        list(Obs=sum(TB.L),
                            Prob=TB.L/sum(TB.L),
                            Pred.class=Pred.class,
                            Dev=Dev.leaf(Lcut[[z]][[1]]),
                            Dev.test=Dev.test,
                            fit.tr=sum(Pred.class == Lcut[[z]][[1]]),
                            fit.test=fit.test)
                    }
                }
            ),
        right=lapply(1:length(Rcut),function(k) {
            if (length(Rcut[[k]][[1]]) > minSplit && levelN > levv)
            {
                sp.sl(Rcut[[k]][[1]],Rcut[[k]],
                if(length(Rtest) > 0) Rtest[[k]],
                        dmin=dmin,
                        minSplit=minSplit,
                        minBucket=minBucket,
                        topn=topn,
                        meth=meth,
                        topnmeth=topnmeth,
                        levelN=levelN,
                        levv=levv,
                        lstep=lstep,
                        tl=tl,
                        kfold=kfold,oldspV=splvar[k])
            } else {
                TB.R <- .Call("tw_table",as.integer(Rcut[[k]][[1]]),levels(Rcut[[k]][[1]]),PACKAGE="TWIX")
                Pred.class <-attr(which.max(TB.R),"names")
                if(length(Rtest) > 0 &&
                    length(Rtest[[k]][[1]]) > 0)
                {
                    Dev.test <- Dev.leaf(Rtest[[k]][[1]])
                    id.test <- Pred.class == Rtest[[k]][[1]]
                    fit.test <- .Call("tw_table",as.integer(Rtest[[k]][[1]][id.test]),
                                levels(Rtest[[k]][[1]][id.test]),PACKAGE="TWIX")
                } else {
                    Dev.test <-fit.test <- 0
                }
                list(Obs=sum(TB.R),
                    Prob=TB.R/sum(TB.R),
                    Pred.class=Pred.class,Dev=Dev.leaf(Rcut[[k]][[1]]),
                    Dev.test=Dev.test,
                    fit.tr=sum(Pred.class == Rcut[[k]][[1]]),
                    fit.test=fit.test)
              }
             }
          )
        )
    }
	else {
            ylev <- .Call("tw_table",as.integer(rsp),levels(rsp),PACKAGE="TWIX")
            Pred.class <-attr(which.max(ylev),"names")
            Dev <- Dev.leaf(rsp)
            fit.tr <- sum(Pred.class == rsp)
            if(length(test[[1]]) > 0) {
                Dev.test <- Dev.leaf(test[[1]])
                id.test <- Pred.class == test[[1]]
                fit.test <- .Call("tw_table",as.integer(test[[1]][id.test]),
                                levels(test[[1]][id.test]),PACKAGE="TWIX")
            } else {
                Dev.test <- fit.test <- 0
                }
            m<-test<-0
            E <- list(Obs=sum(ylev),Prob=ylev/sum(ylev),
                        Pred.class=Pred.class,Dev=Dev,
                        Dev.test=Dev.test,fit.tr=fit.tr,
                        fit.test=fit.test)
		}
    E
    }
    sp.sl(rsp,m,test=test.data,topn=topN)
}

