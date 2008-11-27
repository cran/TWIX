
sp.bagg <- function(rsp,m,test.data=NULL,Devmin=0.01,minsplit=20,
                minbucket=10,cluster=NULL,topN=1,topn.method="complete",
                level=30, method="deviance",tol=0.15, splitf="deviance", st=1, k=0, robust=FALSE)
{
	sp.twix <- function(rsp,m, test.d=NULL, dmin=Devmin, minSplit=minsplit,
                minBucket=minbucket, clname=cluster, topn=topN, topn.meth=topn.method,
                levelN=level,meth=method,tl=tol,splf=splitf,lstep=st,K=k,rob=robust,
				lev=0,cllevel=NULL)
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
			if(robust && lev < 2){
				imp <- impor(m[,2:n],rsp,runs=15)
				TTopn <- 20
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
			S <- lapply(m[2:n],splitt_padj,rsp=as.numeric(rsp))
		}
		if( K != 0 && lev < 2 && nrow(m) > 60) {
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
		else {
			if(S[[1]]$globD > 0){
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
                    test=FALSE,
                    level=lev,K=K)
		if( K != 0 && lev < 2 && nrow(m) > 60) {
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
	
    splvar <- vector()
    Lcut <- Ltest <- vector("list",length(topn))
	Rcut <- Rtest <- vector("list",length(topn))

    make.node <- function(y,dev=S[[1]][y],gdev=S[[2]],
            spoint=S[[3]][[id_Var[y]]],
			var.id=S[[5]][y],id.var=id_Var[y],
            mindev=dmin,minbucket=minBucket,data=m,
			newdata=test.d,k=h,LL=lev) {
        ans <- 0
        if(length(dev) > 0 && dev > 0){
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
                    if(!is.null(newdata))
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
    if(length(S) > 0 && length(S[[1]]) > 0)
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
    S<-id_Var<-0
    splitnode_par <- function(z,Lcut,Ltest,Sval,dmin,minSplit,minBucket,top,
                            meth,topnmeth,levelN,ll,stt,Tol,kfold,splitfun)
    {
        if(length(Lcut[[z]][[1]]) > minSplit && levelN > ll )
        {
            sp.slave(Lcut[[z]][[1]],Lcut[[z]],
                if(!is.null(Ltest[[z]])) Ltest[[z]],
                Dmin=dmin,minsplit=minSplit,
                minbucket=minBucket,topN=top,
                method=meth,topn.method=topnmeth,
                level=levelN,lev=ll,st=stt,tol=Tol,
                K=kfold,splitf=splitfun)
        }
        else {
            TB.L <- table(Lcut[[z]][[1]])
            Pred.class <- names(sort(TB.L,decreasing=TRUE)[1])
            if(!is.null(Ltest[[z]])) {
                Dev.test <- .Call("Dev_leaf",Ltest[[z]][[1]],PACKAGE="TWIX")
                id.test <- Pred.class == Ltest[[z]][[1]]
                fit.test <- table(Ltest[[z]][[1]][id.test])
            }
            else {
                Dev.test <- fit.test <- 0
            }
            list(Obs=sum(TB.L),
                Prob=TB.L/sum(TB.L),
                Pred.class=Pred.class,Dev=.Call("Dev_leaf",Lcut[[z]][[1]],PACKAGE="TWIX"),
                Dev.test=Dev.test,fit.tr=sum(Pred.class ==Lcut[[z]][[1]]),
                fit.test=fit.test)
        }
    }
    if(lev >= cllevel && !is.null(clname)) {
        if (length(Lcut) > 0 && length(Rcut) > 0) {
            m<-test.d<-0
            E <-list(split=Sval,
                left=clusterApplyLB(clname,1:length(Lcut),
                    splitnode_par,Lcut,Ltest,Sval=Sval,minSplit=minSplit,
                    dmin=dmin,minBucket=minBucket,top=topN,meth=meth,
                    topnmeth=topn.meth,levelN=levelN,ll=lev,stt=st,Tol=tol,
                    kfold=K,splitfun=splf),
                right=clusterApplyLB(clname,1:length(Rcut),
                    splitnode_par,Rcut,Rtest,Sval=Sval,minSplit=minSplit,
                    dmin=dmin,minBucket=minBucket,top=topN,meth=meth,
                    topnmeth=topn.meth,levelN=levelN,ll=lev,stt=st,Tol=tol,
                    kfold=K,splitfun=splf)
                )
        }
        else {
            ylev <- table(m[1])
            Pred.class <- names(sort(ylev,decreasing=TRUE)[1])
            Dev <- .Call("Dev_leaf",m[[1]],PACKAGE="TWIX")
            fit.tr <- sum(Pred.class == m[[1]])
            if(length(test.d[[1]]) > 0) {
                Dev.test <- .Call("Dev_leaf",test.d[[1]],PACKAGE="TWIX")
				id.test <- Pred.class == test.d[[1]]
                fit.test <- table(test.d[[1]][id.test])
            }
            else {
                Dev.test <- fit.test <- 0
            }
            m<-test.d<-0
            E <- list(Obs=sum(ylev),Prob=ylev/sum(ylev),
                    Pred.class=Pred.class,Dev=Dev,Dev.test=Dev.test,
                    fit.tr=fit.tr,fit.test=fit.test)
        }
        E
    }
    else {
        if (h > 1) {
            m <- test.d <- 0
            E <-list(split=Sval,
            left=lapply(1:(h-1),
                    function(z) {
                    if(length(Lcut[[z]][[1]]) > minSplit && levelN > lev){
                    sp.twix(Lcut[[z]][[1]],
                        Lcut[[z]],
                        if(length(Ltest) > z) Ltest[[z]],
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
					else {
						yLcut <- Lcut[[z]][[1]]
                        TB.L <- .Call("tw_table",yLcut,PACKAGE="TWIX")
                        Pred.class <- attr(which.max(TB.L),"names")
						if(length(Ltest) > z){
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
        right=lapply(1:(h-1),
            function(k) {
            if(length(Rcut[[k]][[1]]) > minSplit && levelN > lev){
                sp.twix(Rcut[[k]][[1]],
                Rcut[[k]],
                if(length(Rtest) > k) Rtest[[k]],
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
                if(length(Rtest) > k){
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
            if(length(test.d[[1]]) > 0) {
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
        }
    E
  }

	K <- sp.twix(rsp,m)

    t <- 0
    b <- function(n,Svar,Sp,d,ftr,node.cl,L,R){
        tree <- list()
        for (i in 1:length(L)) {
            for (j in 1:length(R)) {
                if(length(L[[i]]$id) > 1 | length(R[[j]]$id) > 1){
                    id <- c(n,L[[i]]$id,R[[j]]$id)
                    Splitvar <- c(Svar,L[[i]]$Splitvar,R[[j]]$Splitvar)
                    Splitp <- c(Sp[1],L[[i]]$Splitp,R[[j]]$Splitp)
                    dev <- d+L[[i]]$dev+R[[j]]$dev
                    fit.tr <- ftr+L[[i]]$fit.tr+R[[j]]$fit.tr
                    tree[[t<-t+1]]<-list(id=id,dev=dev,
                        fit.tr=fit.tr,Splitvar=Splitvar,
                        Splitp=Splitp,Pred.class=node.cl)
                }
#                else if(L[[i]]$Pred.class == R[[j]]$Pred.class) {
#                    id<-0
#                    Splitvar<-Svar
#                    Splitp<-Sp[1]
#                    dev<-0
#                    fit.tr<-ftr+L[[i]]$fit.tr+R[[j]]$fit.tr
#                    tree[[t<-t+1]]<-list(id=id,dev=dev,
#                        fit.tr=fit.tr,Splitvar=Splitvar,
#                        Splitp=Splitp,Pred.class=node.cl)
#                }
                else {
                    id <- c(n,L[[i]]$id,R[[j]]$id)
                    Splitvar <- c(Svar,L[[i]]$Splitvar,R[[j]]$Splitvar)
                    Splitp <- c(Sp[1],L[[i]]$Splitp,R[[j]]$Splitp)
                    dev <- d+L[[i]]$dev+R[[j]]$dev
                    fit.tr <- ftr+L[[i]]$fit.tr+R[[j]]$fit.tr
                    tree[[t<-t+1]]<-list(id=id,dev=dev,
                        fit.tr=fit.tr,Splitvar=Splitvar,
                        Splitp=Splitp,Pred.class=node.cl)
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
                        0,a$split[[k]]$Pred.class,
                    if(length(a$left[[k]]) == 3)
                        s(a$left[[k]])
                    else
                        list(list(id=0,Splitvar=0,Splitp=0,dev=0,
                        fit.tr=a$left[[k]]$fit.tr,
                        Pred.class=a$left[[k]]$Pred.class)),
                    if(length(a$right[[k]]) == 3)
                        s(a$right[[k]])
                    else
                        list(list(id=0,Splitvar=0,Splitp=0,dev=0,
                        fit.tr=a$right[[k]]$fit.tr,
                        Pred.class=a$right[[k]]$Pred.class))
                )
            )
        }
    }
    if(!is.null(K$split)){		
        KK <- s(K)
        id.leng <- sapply(KK,function(x) length(x$id))
        id.out <- which(id.leng == 1)
        if(length(id.out) > 0 && length(id.out) < length(KK)){
            KK <- KK[-id.out]
            id.leng <- id.leng[-id.out]
        }
        if(length(id.out) == length(id.leng)){
            KK<-list(KK[[1]])
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
        l<-list()
        l$dev.test<-l$dev<-l$id<-0
        K<-c(l,K)
        KK<-gr<-list(K)
        l<-list()
        l$split<-list(K)
        K<-l
    }
	if(length(topN) == 1){
		topn <- 1:topN
    }
	else{
		topn <- topN
	}
	if(length(topn) > 1){
		KK <- KK[my.sort(sapply(KK,function(x)x$dev),
			index.return=TRUE,decreasing=TRUE)$ix]
		KK<-KK[my.sort(sapply(KK,function(x)x$fit.tr),
			index.return=TRUE,decreasing=TRUE)$ix]
	}
	
	class(KK) <- "id.tree"
    database <- list(multitree=K,trees=KK)
	if(splitf != "p-adj")
		class(database) <- c("TWIX","deviance")
	else
		class(database) <- c("TWIX","p-adj")
    database
}

