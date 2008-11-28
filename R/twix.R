TWIX <- function(formula,data=NULL,test.data=NULL,subset=NULL,
        method="deviance",topn.method="complete",cluster=NULL,
        minsplit=30,minbucket=round(minsplit/3),Devmin=0.05,
		splitf="deviance",topN=1,level=30,st=1,cl.level=2,
		tol=0.15,score=1,k=0,verbose=FALSE,trace.plot=FALSE,...)
{
    call <- match.call()
    m <- match.call(expand=FALSE)
    m$method <-m$topn.method <- m$cl.level <- m$test.data <- NULL
    m$cluster <- m$minsplit <- m$minbucket <- NULL
    m$Devmin <- m$topN <- m$level <- m$st <- m$tol <- m$... <-NULL
    m$score <- m$k <- m$trace.plot <- m$robust<- m$splitf <-NULL
    m[[1]] <- as.name("model.frame.default")
    m <- eval(m, parent.frame())
	
	data.Na <- apply(m, 1, function(x) any(is.na(x)))
	if(sum(data.Na) > 0){
		stop("\n Missing values in training data! \n\n")
		}
    if(splitf == "p-adj")
        Devmin <- 1-Devmin
    twix.data <- function(Dt){
		fact <- sapply(Dt, function(x) !is.null(levels(x)))
		char <- sapply(Dt, is.character)
		sDt <- 1:ncol(Dt)
		if(any(fact | char)) {
			for (j in sDt[char])
				Dt[,j] <- factor(Dt[[j]])
		}
		Dt
	}
    if(is.null(test.data) || is.data.frame(test.data)){
        m <- twix.data(m)
    }
    else if(is.matrix(test.data)){
        stop("\n   test.data must be a data.frame!\n")
        }
    else if(length(test.data) == 1){
        test.data <- NULL
        }
    if(!is.null(test.data)) {
        test.data <- model.frame(formula,test.data)
        ntest <- dim(test.data)[1]
    }
	else{
        ntest <- 1
	}
    ntr <- dim(m)[1]
    rsp <- model.extract(m, "response")
    robust <- bag <- FALSE
    j <- i <- 1
    if(!is.null(names(list(...)))){
        if(names(list(...)) == "bag")
			bag <- TRUE
		if(names(list(...)) == "robust")
			robust <- TRUE
	}
    if(method=="grid" && topn.method != "single")
        topn.method <- "single"
    if(!is.null(cluster)) {
        clusterSetupRNG(cluster)
        clusterEvalQ(cluster, library(TWIX))
        }
    sp <- function(rsp,m,test.d=test.data,dmin=Devmin,minSplit=minsplit,
                minBucket=minbucket,clname=cluster,cllevel=cl.level,topn=topN,topn.meth=topn.method,
                levelN=level,lev=0,meth=method,lstep=st,tl=tol,K=k,rob=robust,splf=splitf)
    {
    n <- length(m)
    E <- vector("list",3)
	Sval <-list()
    lev <- lev+1
    if (length(topN) == 1 ) {
        topn <- 1:topN
    }
    else {
        if (!is.na(topN[lev]) && topN[lev] != 1) {
            topn <- 1:topN[lev]
        }
        else
            topn <- 1
    }
    if(lev == 2 && meth == "grid") {
        meth <- "deviance"
        topn.meth <- "complete"
        }
    if(lev > 2 && meth == "local"){
        meth <- "deviance"
        }
    k.topn <- round(length(topn)/(n-1))
    if(k.topn == 0) k.topn <-1
    if (is.null(clname)) {
		if(splf=="deviance"){
			if(robust && lev < 3){
				imp<-impor(m[,2:n],rsp,runs=10)
				TTopn<-20
				S <- lapply(m[2:n],
                            splitt,
                            rsp,
                            meth=meth,
                            lstep=lstep,
                            topn=TTopn
                            ,topn.meth="complete",
                            test=FALSE,
                            level=lev,K=K)
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
                            level=lev,K=K)
			}
		}
		if(splf == "p-adj"){
			S <- lapply(m[2:n],
				splitt_padj,
				rsp=as.numeric(rsp))
		}
		if(splf != "p-adj" & splf != "deviance")
			stop("\n   splitf must be deviance or p-adj! \n")
		if( K != 0 && lev < 3 && nrow(m) > 60) {
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
    }
	else {
		if(splf == "deviance"){
				if(length(topn) == 1)
					Ttest <- TRUE
				else
					Ttest <- FALSE
			S <- clusterApplyLB(clname,
                    m[2:n],
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
                    level=lev,K=K)
		}
		if(splf == "p-adj"){
			S <- clusterApplyLB(clname,
                    m[2:n],
                    splitt_padj,
                    rsp=as.numeric(rsp))
		}
		if(splf != "p-adj" & splf != "deviance")
			stop("\n   splitf must be deviance or p-adj!\n")
		if( K != 0 && lev < 3 && nrow(m) > 60) {
			if (S[[1]]$globD != 0){
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
		else {
			if (S[[1]]$globD != 0){
				k <- 0
				if (topn.meth != "single")  k <- 1
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
    }
    splvar <- vector()
    Lcut <- Ltest <- vector("list",length(topn))
	Rcut <- Rtest <- vector("list",length(topn))
	
	make.node <- function(y,dev=S[[1]][y],gdev=S[[2]],
            spoint=S[[3]][[id_Var[y]]],var.id=S[[5]][y],
			id.var=id_Var[y],mindev=dmin,minbucket=minBucket,
			data=m,newdata=test.d,k=h,LL=lev) {
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
    splitnode_par <- function(z,Lcut,Ltest,Sval,dmiN,minSplit,minBucket,top,
                            meth,topnmeth,levelN,ll,stt,Tol,kfold,splitfun)
    {
        if(length(Lcut[[z]][[1]]) >= minSplit && levelN > ll )
        {
            sp.slave(Lcut[[z]][[1]],Lcut[[z]],
                if(!is.null(Ltest[[z]])) Ltest[[z]],
                Dmin=dmiN,minsplit=minSplit,
                minbucket=minBucket,topN=top,
                method=meth,topn.method=topnmeth,
                level=levelN,lev=ll,st=stt,tol=Tol,
                K=kfold,splitf=splitfun)
        }
        else {
			yLcut <- Lcut[[z]][[1]]
            TB.L <- .Call("tw_table",yLcut,PACKAGE="TWIX")
            Pred.class <- names(sort(TB.L,decreasing=TRUE)[1])
            if(length(Ltest) >= z && !is.null(Ltest[[z]][[1]])){
				yLtest <- Ltest[[z]][[1]]
                Dev.test <- .Call("Dev_leaf",yLtest,PACKAGE="TWIX")
                fit.test <- sum(Pred.class == yLtest)
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
    if(lev >= cllevel && !is.null(clname)) {
        if (h > 1) {
            E <-list(split=Sval,
                left=clusterApplyLB(clname,1:(h-1),
                    splitnode_par,Lcut,Ltest,Sval=Sval,dmiN=dmin,minSplit=minSplit,
                    minBucket=minBucket,top=topN,meth=meth,
                    topnmeth=topn.meth,levelN=levelN,ll=lev,stt=st,Tol=tol,
                    kfold=K,splitfun=splf),
                right=clusterApplyLB(clname,1:(h-1),
                    splitnode_par,Rcut,Rtest,Sval=Sval,dmiN=dmin,minSplit=minSplit,
                    minBucket=minBucket,top=topN,meth=meth,
                    topnmeth=topn.meth,levelN=levelN,ll=lev,stt=st,Tol=tol,
                    kfold=K,splitfun=splf)
                )
        }
        else {
            ylev <- .Call("tw_table",rsp,PACKAGE="TWIX")
            Pred.class <- attr(which.max(ylev),"names")
            Dev <- .Call("Dev_leaf",rsp,PACKAGE="TWIX")
            fit.tr <- sum(Pred.class == rsp)
            if(length(test.d[[1]]) > 0) {
				ytest <- test.d[[1]]
                Dev.test <- .Call("Dev_leaf",ytest,PACKAGE="TWIX")
                fit.test <- sum(Pred.class == ytest)
            }
            else {
                Dev.test <- fit.test <- 0
            }
            E <- list(Obs=sum(ylev),
					Prob=ylev/sum(ylev),
                    Pred.class=Pred.class,
					Dev=Dev,Dev.test=Dev.test,
                    fit.tr=fit.tr,
					fit.test=fit.test)
        }
        E
    }
    else {
        if(h > 1){
            E <- list(split=Sval,
            left=lapply(1:(h-1), function(z){
                    if(length(Lcut[[z]][[1]]) >= minSplit && levelN > lev){
                    sp(Lcut[[z]][[1]],
                        Lcut[[z]],
                        if(length(Ltest) >= z && !is.null(Ltest[[z]])) Ltest[[z]],
                        dmin=dmin,
                        minSplit=minSplit,
                        minBucket=minBucket,
                        clname=clname,
                        cllevel=cllevel,
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
						#rm(list=c("m","test.d"))
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
			right=lapply(1:(h-1),
				function(k) {
				if(length(Rcut[[k]][[1]]) >= minSplit && levelN > lev )
				{
                sp(Rcut[[k]][[1]],
					Rcut[[k]],
					if(length(Rtest) >= k && !is.null(Rtest[[k]])) Rtest[[k]],
					dmin=dmin,
					minSplit=minSplit,
					minBucket=minBucket,
					clname=clname,
					cllevel=cllevel,
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
                    Pred.class=Pred.class,Dev=.Call("Dev_leaf",yRcut,PACKAGE="TWIX"),
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
                id.test <- Pred.class == ytest
                fit.test <- sum(id.test)
				rm(list="test.d")
			}
			else{
                Dev.test <- fit.test <- 0
			}
            rm(list=c("m"))
            E <- list(Obs=sum(ylev),Prob=ylev/sum(ylev),
                        Pred.class=Pred.class,Dev=Dev,
                        Dev.test=Dev.test,fit.tr=fit.tr,
                        fit.test=fit.test)
		}
	}
    E
  }
    K <- sp(rsp,m,test.data)
    t <- 0
    b <- function(n,Svar,Sp,d,d.t,ftr,ftest,node.cl,sdtr,Obsn,dis,kst,ddD,L,R){
        tree <- list()
        for (i in 1:length(L)) {
            for (j in 1:length(R)) {
                if(length(L[[i]]$id) > 1 | length(R[[j]]$id) > 1){
                    id <- c(n,L[[i]]$id,R[[j]]$id)
                    Splitvar <- c(Svar,L[[i]]$Splitvar,R[[j]]$Splitvar)
                    Splitp <- c(Sp[1],L[[i]]$Splitp,R[[j]]$Splitp)
                    dev <- d+L[[i]]$dev+R[[j]]$dev
                    fit.tr <- ftr+L[[i]]$fit.tr+R[[j]]$fit.tr
                    fit.test <- ftest+L[[i]]$fit.test+R[[j]]$fit.test
                    dev.test<- L[[i]]$dev.test+R[[j]]$dev.test
                    sd.tr <- c(sdtr,L[[i]]$sd.tr,R[[j]]$sd.tr)
                    Obs <- c(Obsn,L[[i]]$Obs,R[[j]]$Obs)
                    dist <- c(dis,L[[i]]$dist,R[[j]]$dist)
                    ks.t <- c(kst,L[[i]]$ks.t,R[[j]]$ks.t)
					dD <- c(ddD,L[[i]]$dD,R[[j]]$dD)
                    tree[[t<-t+1]]<-list(id=id,dev=dev,dev.test=dev.test,
                        fit.tr=fit.tr,fit.test=fit.test,Splitvar=Splitvar,
                        Splitp=Splitp,Pred.class=node.cl,sd.tr=sd.tr,
						Obs=Obs,dist=dist,ks.t=ks.t,dD=dD)
                }
                else if(L[[i]]$Pred.class == R[[j]]$Pred.class) {
                    id<-0
                    Splitvar<-Svar
                    Splitp<-Sp[1]
                    dev<-0
                    fit.tr<-ftr+L[[i]]$fit.tr+R[[j]]$fit.tr
                    fit.test<-ftest+L[[i]]$fit.test+R[[j]]$fit.test
                    dev.test<-d.t
                    sd.tr <- sdtr
                    Obs <- L[[i]]$Obs+R[[j]]$Obs
                    dist <- ks.t<-0
                    tree[[t<-t+1]]<-list(id=id,dev=dev,dev.test=dev.test,
                        fit.tr=fit.tr,fit.test=fit.test,Splitvar=Splitvar,
                        Splitp=Splitp,Pred.class=node.cl,sd.tr=sd.tr,
						Obs=Obs,dist=dist,ks.t=ks.t,dD=0)
                }
                else {
                    id <- c(n,L[[i]]$id,R[[j]]$id)
                    Splitvar <- c(Svar,L[[i]]$Splitvar,R[[j]]$Splitvar)
                    Splitp <- c(Sp[1],L[[i]]$Splitp,R[[j]]$Splitp)
                    dev <- d+L[[i]]$dev+R[[j]]$dev
                    fit.tr <- ftr+L[[i]]$fit.tr+R[[j]]$fit.tr
                    fit.test <- ftest+L[[i]]$fit.test+R[[j]]$fit.test
                    dev.test<- L[[i]]$dev.test+R[[j]]$dev.test
                    sd.tr <- c(sdtr,L[[i]]$sd.tr,R[[j]]$sd.tr)
                    Obs <- c(Obsn,L[[i]]$Obs,R[[j]]$Obs)
                    dist <- c(dis,L[[i]]$dist,R[[j]]$dist)
                    ks.t <- c(kst,L[[i]]$ks.t,R[[j]]$ks.t)
					dD <- c(ddD,L[[i]]$dD,R[[j]]$dD)
                    tree[[t<-t+1]]<-list(id=id,dev=dev,dev.test=dev.test,
                        fit.tr=fit.tr,fit.test=fit.test,Splitvar=Splitvar,
                        Splitp=Splitp,Pred.class=node.cl,sd.tr=sd.tr,
						Obs=Obs,dist=dist,ks.t=ks.t,dD=dD)
                }
            }
        }
    tree
    }
    s <- function(a) {
        id <- list()
        for (k in 1:length(a$split)) {
            id <- c(id,
                b(k,a$split[[k]]$Splitvar,a$split[[k]]$Splitp,a$split[[k]]$Dev,
                        a$split[[k]]$Dev.test,0,0,
                        a$split[[k]]$Pred.class,
                        max(a$split[[k]]$Prob),0,a$split[[k]]$dist,
                        a$split[[k]]$ks.t,
						a$split[[k]]$dD,
                    if(length(a$left[[k]]) == 3)
                        s(a$left[[k]])
                    else
                        list(list(id=0,Splitvar=0,Splitp=0,dev=0,
                        dev.test=a$left[[k]]$Dev.test,
                        fit.tr=a$left[[k]]$fit.tr,
                        fit.test=a$left[[k]]$fit.test,
                        Pred.class=a$left[[k]]$Pred.class,
                        sd.tr=max(a$left[[k]]$Prob),
						Obs=a$left[[k]]$Obs,
						dist=0,ks.t=0,dD=0)),
                    if(length(a$right[[k]]) == 3)
                        s(a$right[[k]])
                    else
                        list(list(id=0,Splitvar=0,Splitp=0,dev=0,
                        dev.test=a$right[[k]]$Dev.test,
                        fit.tr=a$right[[k]]$fit.tr,
                        fit.test=a$right[[k]]$fit.test,
                        Pred.class=a$right[[k]]$Pred.class,
                        sd.tr=max(a$right[[k]]$Prob),
						Obs=a$right[[k]]$Obs,
						dist=0,ks.t=0,dD=0))
                )
            )
        }
    }
	
    b_padj <- function(n,Svar,Sp,d,d.t,ftr,ftest,node.cl,sdtr,Obsn,dis,kst,ddD,L,R){
        tree <- list()
        for (i in 1:length(L)) {
            for (j in 1:length(R)) {
                if(length(L[[i]]$id) > 1 | length(R[[j]]$id) > 1){
                    id <- c(n,L[[i]]$id,R[[j]]$id)
                    Splitvar <- c(Svar,L[[i]]$Splitvar,R[[j]]$Splitvar)
                    Splitp <- c(Sp[1],L[[i]]$Splitp,R[[j]]$Splitp)
                    dev <- L[[i]]$dev+R[[j]]$dev
                    fit.tr <- ftr+L[[i]]$fit.tr+R[[j]]$fit.tr
                    fit.test <- ftest+L[[i]]$fit.test+R[[j]]$fit.test
                    dev.test<- L[[i]]$dev.test+R[[j]]$dev.test
                    sd.tr <- c(sdtr,L[[i]]$sd.tr,R[[j]]$sd.tr)
                    Obs <- c(Obsn,L[[i]]$Obs,R[[j]]$Obs)
                    dist <- c(dis,L[[i]]$dist,R[[j]]$dist)
                    ks.t <- c(kst,L[[i]]$ks.t,R[[j]]$ks.t)
					dD <- c(ddD,L[[i]]$dD,R[[j]]$dD)             
                    tree[[t<-t+1]]<-list(id=id,dev=dev,dev.test=dev.test,
                        fit.tr=fit.tr,fit.test=fit.test,Splitvar=Splitvar,
                        Splitp=Splitp,Pred.class=node.cl,sd.tr=sd.tr,Obs=Obs,dist=dist,ks.t=ks.t,dD=dD)
                }
                else if(L[[i]]$Pred.class == R[[j]]$Pred.class) {
                    id<-0
                    Splitvar<-Svar
                    Splitp<-Sp[1]
                    dev<-d
                    fit.tr<-ftr+L[[i]]$fit.tr+R[[j]]$fit.tr
                    fit.test<-ftest+L[[i]]$fit.test+R[[j]]$fit.test
                    dev.test<-d.t
                    sd.tr <- sdtr
                    Obs <- L[[i]]$Obs+R[[j]]$Obs
                    dist <- ks.t<-0
                    tree[[t<-t+1]]<-list(id=id,dev=dev,dev.test=dev.test,
                        fit.tr=fit.tr,fit.test=fit.test,Splitvar=Splitvar,
                        Splitp=Splitp,Pred.class=node.cl,sd.tr=sd.tr,Obs=Obs,dist=dist,ks.t=ks.t,dD=0)
                }
                else {
                    id <- c(n,L[[i]]$id,R[[j]]$id)
                    Splitvar <- c(Svar,L[[i]]$Splitvar,R[[j]]$Splitvar)
                    Splitp <- c(Sp[1],L[[i]]$Splitp,R[[j]]$Splitp)
                    dev <- L[[i]]$dev+R[[j]]$dev
                    fit.tr <- ftr+L[[i]]$fit.tr+R[[j]]$fit.tr
                    fit.test <- ftest+L[[i]]$fit.test+R[[j]]$fit.test
                    dev.test<- L[[i]]$dev.test+R[[j]]$dev.test
                    sd.tr <- c(sdtr,L[[i]]$sd.tr,R[[j]]$sd.tr)
                    Obs <- c(Obsn,L[[i]]$Obs,R[[j]]$Obs)
                    dist <- c(dis,L[[i]]$dist,R[[j]]$dist)
                    ks.t <- c(kst,L[[i]]$ks.t,R[[j]]$ks.t)
					dD <- c(ddD,L[[i]]$dD,R[[j]]$dD)                    
                    tree[[t<-t+1]]<-list(id=id,dev=dev,dev.test=dev.test,
                        fit.tr=fit.tr,fit.test=fit.test,Splitvar=Splitvar,
                        Splitp=Splitp,Pred.class=node.cl,sd.tr=sd.tr,Obs=Obs,dist=dist,ks.t=ks.t,dD=dD)
                }
            }
        }
    tree
    }

	s_padj <- function(a) {
        id <- list()
        for (k in 1:length(a$split)) {
            id <- c(id,
                b_padj(k,a$split[[k]]$Splitvar,a$split[[k]]$Splitp,a$split[[k]]$dist,
                        a$split[[k]]$Dev.test,0,0,
                        a$split[[k]]$Pred.class,
                        max(a$split[[k]]$Prob),0,a$split[[k]]$dist,
                        a$split[[k]]$ks.t,
						a$split[[k]]$dD,
                    if(length(a$left[[k]]) == 3)
                        s_padj(a$left[[k]])
                    else
                        list(list(id=0,Splitvar=0,Splitp=0,dev=a$left[[k]]$Dev,
                        dev.test=a$left[[k]]$Dev.test,
                        fit.tr=a$left[[k]]$fit.tr,
                        fit.test=a$left[[k]]$fit.test,
                        Pred.class=a$left[[k]]$Pred.class,
                        sd.tr=max(a$left[[k]]$Prob),
						Obs=a$left[[k]]$Obs,
						dist=0,ks.t=0,dD=0)),
                    if(length(a$right[[k]]) == 3)
                        s_padj(a$right[[k]])
                    else
                        list(list(id=0,Splitvar=0,Splitp=0,dev=a$right[[k]]$Dev,
                        dev.test=a$right[[k]]$Dev.test,
                        fit.tr=a$right[[k]]$fit.tr,
                        fit.test=a$right[[k]]$fit.test,
                        Pred.class=a$right[[k]]$Pred.class,
                        sd.tr=max(a$right[[k]]$Prob),
						Obs=a$right[[k]]$Obs,
						dist=0,ks.t=0,dD=0))
                )
            )
        }
    }
    DF <- id.agr <- agr.id <- 1
    if(!is.null(K$split)){
		if(splitf != "p-adj"){
			KK <- s(K)
		}
		else{
			KK <- s_padj(K)
		}
        id.leng <- sapply(KK,function(x)length(x$id))
        id.out <- which(id.leng == 1)
        no.score <- FALSE
        if(length(id.out) > 0 && length(id.out) < length(KK)){
            KK <- KK[-id.out]
            id.leng <- id.leng[-id.out]
        }
        if(length(id.out) == length(id.leng)){
            KK <- list(KK[[1]])
            gr <- list(KK[[1]])
            no.score <- TRUE
            gr.id <- tw.tic <- gr.tic <- 1
        }
        v.fit <- function(x,y){
            y-x$fit.test
        }
        tr.fit <- function(x,z){
            z-x$fit.tr
        }
        for(i in 1:length(KK)) {
            KK[[i]]$dev.test <- K$split[[1]]$Dev.test-KK[[i]]$dev.test
            KK[[i]]$fit.tr <- sum(KK[[i]]$fit.tr)/ntr
            KK[[i]]$fit.test <- sum(KK[[i]]$fit.test)/ntest
        }
    }
    else{
		l <- list(dev.test=0, dev=0, id=0)
        K <- c(l,K)
        KK <- gr <- list(K)
        l <- list()
        l$split <- list(K)
        K <- l
        no.score <- TRUE
        gr.id <- tw.tic <- gr.tic <- agr.id <- 1
    }
    if(!no.score){
		gr <- list(KK[[1]])
		gr.l <- length(gr[[1]]$id[gr[[1]]$id==0])
		if(method == "grid") {
			gr <- list(s(sp(rsp,m,test.data,topn=1,meth="deviance"))[[1]])
			gr[[1]]$fit.tr <- gr[[1]]$fit.tr/ntr
			gr[[1]]$fit.test <- gr[[1]]$fit.test/ntest
		}
		cht <- function(x,n=1,levell=3){
			ansl<-colour<-fit<-out<-vector()
			y<-x
			x<-x$id
			l<-r<-j<-0
			left<-TRUE
			k<-1
			for(i in 1:length(x)){
				if(i < levell+2 && x[i] != 0 && left){
					ansl[k]<-y$Splitvar[i]
					out[k]<-y$Splitp[i]
					fit[k]<-y$fit.test
					colour[k]<-k
					k<-k+1
				}
				if(x[i] == 0 && left){     
					l<-l+1
					if(i-l == l) left<-FALSE
				}
				else if(x[i]==0){
					r<-r+1
				}
				if(j < levell-1 && x[i] != 0 && !left){
					j<-j+1
					ansl[k]<-y$Splitvar[i]
					out[k]<-y$Splitp[i]
					fit[k]<-y$fit.test
					colour[k]<-j+1
					k<-k+1
				}
			}
			data.frame(ansl,out,colour,fit)
		}
		max.fit <- gr.tic <- tw.tic <- DF <- 0
		form <- function(obj){
			Var<-splitp<-fit<-colour<-vector()
			for(i in 1:length(obj)){
				Var <- c(Var,as.character(obj[[i]][,1]))
				splitp <- c(splitp,obj[[i]][,2])
				colour <- c(colour,obj[[i]][,3])
				fit <- c(fit,obj[[i]][,4])
			}
			fit <- unclass(as.factor(fit))
			fit <- max(fit)-(fit-1)
			max.fit <<- max(fit)
			data <- data.frame(splitp,colour,fit)
			split(data,as.factor(Var))
		}
		if(trace.plot){
			bal <- lapply(KK,function(x) cht(x,levell=2))
			btr <- form(bal)
			outpl <- list(m[-1],btr,max.fit)
			class(outpl) <- "TWIX"
			Devplot(rsp,outpl)
		}
    if(score == 1){
        fit.tr_certaintly <- sapply(KK,function(x) x$fit.tr*(prod(x$sd.tr[x$Obs!=0])))
        gr.tic <- fit.tr_certaintly[1] 
        klm <- my.sort(fit.tr_certaintly,index.return=TRUE,decreasing=TRUE)$ix
        tw.tic <- fit.tr_certaintly[klm[1]]
        KK <- KK[klm]
        fit.tr_certaintly <- klm <- 0
    }
    if(score == 2){
        fit.tr <- sapply(KK,function(x) x$fit.tr)
        dev.tr <- sapply(KK,function(x) x$dev)
        sd.tr <- lapply(KK,function(x) x$sd.tr[x$Obs != 0])
        Nbuck <- sapply(sd.tr,function(x) length(x))
        sd.tr1 <- sapply(sd.tr,function(x) median(x))
        sd.tr <- sapply(sd.tr,function(x) IQR(x))
        sd.tr <- 0.6*sd.tr1/(1+max(sd.tr1))+0.3*(1-sd.tr/(1+max(sd.tr)))
        sd.tr <- sd.tr/max(sd.tr)
        Obs <- sapply(KK,function(x) min(x$Obs[x$Obs != 0]))
        d <- c(which(fit.tr < quantile(fit.tr)[2]),which(dev.tr < quantile(dev.tr)[2]))
        a <- boxplot(Nbuck,plot=FALSE)
        if(minbucket >= 7){
            if(length(KK)> 30)
                id.out <- c(which(Nbuck > a$stats[5,]),which(Nbuck <= a$stats[1,]))
            else
                id.out<-rep(FALSE,length(KK))
            id.agr<-0.6*sd.tr*(1-Nbuck/(1+max(Nbuck)))+0.2*fit.tr/max(fit.tr)+0.2*(1-Obs/max(Obs))
            id.agr[id.out]<-0.1*(id.agr[id.out])
            agr.id <-my.sort(id.agr,index.return=TRUE,decreasing=TRUE)$ix
            Wbuk<-sd.tr*(1-Nbuck/(1+max(Nbuck)))
            Wbuk<-Wbuk/max(Wbuk)
            if(length(KK)> 30)
                Wbuk[id.out]<-(Wbuk[id.out])*0.25
            DF<-0.15*Wbuk+0.6*fit.tr/max(fit.tr)+0.25*(dev.tr/max(dev.tr))
        }
        if(minbucket < 7){
            if(length(KK)> 30)
                id.out <- c(which(Nbuck <= a$stats[1,]))
            else
                id.out<-rep(FALSE,length(KK))
            id.agr<-0.6*sd.tr*(1-Nbuck/(1+max(Nbuck)))+0.2*fit.tr/max(fit.tr)+0.2*(1-Obs/max(Obs))
            id.agr[id.out]<-0.1*(id.agr[id.out])
            agr.id <-my.sort(id.agr,index.return=TRUE,decreasing=TRUE)$ix
            Wbuk<-sd.tr*(1-Nbuck/(1+max(Nbuck)))
            Wbuk<-Wbuk/max(Wbuk)
            if(length(KK)> 30)
                Wbuk[id.out]<-(Wbuk[id.out])*0.25
            DF<-0.1*Wbuk+0.7*fit.tr/max(fit.tr)+0.2*(dev.tr/max(dev.tr))
            }
        gr.tic <- DF[1] 
        if(length(KK)> 20)
            DF[d]<-0.1*(DF[d])
        klm<-my.sort(DF,index.return=TRUE,decreasing=TRUE)$ix
        tw.tic <- DF[klm[1]]
        KK <- KK[klm]
        fit.tr<-dev.tr<-klm<-a<-b<-d<-0
    }
    if(score == 3){
        gr.tic <- KK[[1]]$fit.test
        if(is.data.frame(test.data)){
            KK <- KK[my.sort(sapply(KK,function(x)x$dev),
                index.return=TRUE,decreasing=TRUE)$ix]
            KK<-KK[my.sort(sapply(KK,function(x)x$fit.test),
                index.return=TRUE,decreasing=TRUE)$ix]
        }
        else {
            KK <- KK[my.sort(sapply(KK,function(x)x$dev),
                index.return=TRUE,decreasing=TRUE)$ix]
            KK <- KK[my.sort(sapply(KK,function(x)x$dev.test),
                index.return=TRUE,decreasing=TRUE)$ix]
            KK<-KK[my.sort(sapply(KK,function(x)x$fit.test),
                index.return=TRUE,decreasing=TRUE)$ix]
            }
        id.agr<-1:length(KK)
        tw.tic <- KK[[1]]$fit.test
    }
    if(score == 4){
        gr.tic <- KK[[1]]$fit.test
        id.agr <- 1:length(KK)
        tw.tic <- KK[[1]]$fit.test
	}
    gr.id <- which(sapply(KK,function(x)x$dev) == gr[[1]]$dev)[1]
    }
    if(!verbose || !bag) {
        cat("n =",length(KK),"\n")
        cat("Deviance gain and TIC of the best TWIX-tree:","    ",c(KK[[1]]$dev,tw.tic),"\n")
        if(!is.na(gr.id))
            cat(paste("Deviance gain and TIC of the greedy tree(Nr.",gr.id,"):"," ",sep =""),
                c(gr[[1]]$dev,gr.tic),"\n")
        else
            cat(paste("Deviance gain and TIC of the greedy tree:  "," ",sep = ""),
                c(gr[[1]]$dev,gr.tic),"\n")
    }
    greedy.tree <- list(greedy.tree=gr,id=if(!is.na(gr.id)) gr.id else 0)
    d.range <- sapply(m[-1],function(x)if(is.character(x) | is.factor(x)) length(table(x)) else max(x))
    class(KK) <- class(greedy.tree) <- "id.tree"
    database <- list(formula=formula(terms(formula,data=m)),
        call=call, greedy.tree=greedy.tree, max.range=d.range,
        multitree=K, trees=KK, agg.id=agr.id,
		dist=dist, score=DF)
	if(splitf != "p-adj")
		class(database) <- c("TWIX","deviance")
	else
		class(database) <- c("TWIX","p-adj")
    database
}
