TWIX <- function(formula,data=NULL,test.data=0,subset=NULL,
        method="deviance",topn.method="complete",cluster=NULL,
        minsplit=30,minbucket=round(minsplit/3),Devmin=0.05,
        topN=1,level=30,st=1,cl.level=2,tol=0.15,score=1,k=0,trace.plot=FALSE,...)
{
    call <-match.call()
    m <- match.call(expand=FALSE)
    m$method <-m$topn.method <- m$cl.level <- m$test.data <- NULL
    m$x <- m$y <- m$cluster <- m$minsplit <- m$minbucket <- NULL
    m$Devmin <- m$topN <- m$level <- m$st <- m$tol <- m$... <-NULL
    m$score <- m$k <- m$trace.plot <- NULL
    m[[1]] <- as.name("model.frame.default")
    m <- eval(m, parent.frame())
    icv<-0
    twix.data <- function(Dt,test.data){
        Terms <- attr(Dt,"terms")
        if(is.null(test.data)) {            
            n<-nrow(Dt)
            fold <- 10
            minbuck <-minbucket
            runs<-length(minbuck)
            CV <- array(0,c(fold,runs,5))
            xgr <- 1:fold
            idd<-vector()
            fmla <- as.formula(paste(names(Dt[1])," ~ ", paste(names(Dt[-1]), collapse= "+")))
            nm <- nrow(table(Dt[,1]))
            for(l in 1:5){
                ids <- sample(rep(xgr, length = n), n)
                for(i in 1:runs) {
                    for(j in xgr) { 
                        test <- ids == j
                        test.id <- which(ids == j)
                        train.id <- which(ids == i)
                        train <- !test 
                        rpt <- rpart(fmla,data=Dt[train,]
                            ,parms=list(split="information"),
                            control=rpart.control(cp=0.01,
                            minbucket=minbuck[i],minsplit=3*minbuck[i]))
                        idd<-c(idd,test.id[(predict(rpt,newdata=Dt[test,], type="class")!= Dt[test,1])])
                        #idd<-c(idd,train.id[(predict(rpt,newdata=Dt[train,], type="class")!= Dt[train,1])])
                        conf <- as.matrix(table(predict(rpt,newdata=Dt[test,], type="class"),Dt[test,1])) 
                        CV[j,i,l] <- sum(conf - diag(diag(conf), nm,nm)) 
                    } 
                }
            }
            m1<-apply(CV,2,median)
            minsplit<<-3*minbucket
            CV<-rpt<-test.id<-test<-conf<-0
            idd2<-as.numeric(names(table(idd)))
            icv<<-idd2[table(idd) >= 5]
            cat(paste("Bad observations:  n =",length(icv),"\n"))
            ict<-setdiff(1:nrow(Dt),icv)
            fact <- sapply(Dt, function(x) !is.null(levels(x)))
            char <- sapply(Dt, is.character)
            sDt <- 1:ncol(Dt)
            if(any(fact | char)) {
                for (j in sDt[char])
                    Dt[,j] <- as.factor(Dt[[j]])
            fact <- fact | char
            cat.levels <- lapply(Dt[fact], levels)
            test.data <- Dt[icv,]
            attr(test.data,"levels") <- cat.levels
            Dt <- Dt[ict,]
            attr(Dt, "levels") <- cat.levels
            attr(Dt, "terms") <- Terms
            list(Dt,test.data)
            }
        }
        else{
            fact <- sapply(Dt, function(x) !is.null(levels(x)))
            char <- sapply(Dt, is.character)
            sDt <- 1:ncol(Dt)
            if(any(fact | char)) {
                for (j in sDt[char])
                    Dt[,j] <- as.factor(Dt[[j]])
            fact <- fact | char
            cat.levels <- lapply(Dt[fact], levels)
            attr(test.data,"levels") <- cat.levels
            attr(test.data, "terms") <- Terms
            attr(Dt,"levels") <- cat.levels
            attr(Dt, "terms") <- Terms
            list(Dt,test.data)
            }
        }
    }
    if(is.null(test.data) || is.data.frame(test.data)){
        M <- twix.data(m,test.data)
        m <- M[[1]]
        test.data <- M[[2]]
        rm(list="M")
    }
    else if(is.matrix(test.data)){
        stop("   test.data must be a data.frame!")
        }
    else if(length(test.data) == 1){
        test.data<-NULL
        }
    if(!is.null(test.data)) {
        test.data<-model.frame(formula,test.data)
        ntest <- dim(test.data)[1]
    } else {
        ntest <- 1
        }
    ntr <- dim(m)[1]
    rsp <- model.extract(m, "response")
    bag <- FALSE
    j <- i <- 1
    if(!is.null(names(list(...))) &&
        names(list(...)) == "bag" ) bag <- TRUE
    if(method=="grid" && topn.method != "single")
        topn.method <- "single"
    if(!is.null(cluster)) {
        clusterSetupSPRNG(cluster)
        clusterEvalQ(cluster, library(TWIX))
        }
    sp <- function(rsp,m,test.d=test.data,dmin=Devmin,minSplit=minsplit,
                minBucket=minbucket,clname=cluster,cllevel=cl.level,topn=topN,topn.meth=topn.method,
                levelN=level,lev=0,meth=method,lstep=st,tl=tol,K=k,oldspvar=0)
    {
    n <- length(m)
    E <- Sval <- list()
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
    k.topn <- round(length(topn)/(n-1))
    if(k.topn == 0) k.topn <-1
    Dev.leaf <- function(x) {
        TD <- 0
        CCR <- .Call("tw_table",as.integer(x),levels(x),PACKAGE="TWIX")
        s <- sum(CCR)
        TD <- -sum(s*sapply(CCR,function(x,y){ if(x/y == 0) 0 else (x/y)*log(x/y)},y=s))
        round(TD,digits=6)
    }
    if (is.null(clname)) {
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
                    test=FALSE,
                    level=lev,K=K)
    if( K != 0 && lev < 2 && nrow(m) > 60) {
        S_summ <- Var_id <- id_Var <- vector()
        which <- list()
        globD <-k<-0
        if (S[[1]]$globD != 0){
            if (topn.meth !="single")  k <- 1
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
            if (topn.meth !="single")  k <- 1
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
    } else {
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
        S_summ <- Var_id <- id_Var <- vector()
        which <- list()
        globD <- k <- 0
        if (S[[1]]$globD != 0){
            if (topn.meth !="single")  k <- 1
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
            if (topn.meth !="single")  k <- 1
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
    }
    splvar <- vector()
    Lcut <- Ltest <- Rcut <- Rtest <- list()
    make.node<- function(y,dev=S_summ[[1]][y],gdev=S_summ[[2]],
            spoint=which[[id_Var[y]]],var.id=Var_id[y],id.var=id_Var[y],
            mindev=dmin,minbucket=minBucket,data=m,newdata=test.d,k=h,oldspv=oldspvar,LL=lev) {
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
    S_summ<-which<-Var_id<-id_Var<-0
    splitnode_par <- function(z,Lcut,Ltest,Sval,dmin,minSplit,minBucket,top,
                            meth,topnmeth,levelN,ll,stt,Tol,kfold,oldspV)
    {
        if(length(Lcut[[z]][[1]]) > minSplit && levelN > ll )
        {
            sp.slave(Lcut[[z]][[1]],Lcut[[z]],
                if(length(Ltest) > 0) Ltest[[z]],
                Dmin=dmin,minsplit=minSplit,
                minbucket=minBucket,topN=top,
                method=meth,topn.method=topnmeth,
                level=levelN,lev=ll,st=stt,tol=Tol,
                K=kfold,oldspvar=oldspV[z])
        }
        else {
            TB.L <- table(Lcut[[z]][[1]])
            Pred.class<-names(sort(TB.L,
                decreasing=TRUE)[1])
            if(length(Ltest) > 0 && length(Ltest[[z]][[1]]) > 0) {
                Dev.test <- Dev.leaf(Ltest[[z]][[1]])
                id.test <- Pred.class == Ltest[[z]][[1]]
                fit.test <- table(Ltest[[z]][[1]][id.test])
            }
            else {
                Dev.test <- fit.test <- 0
            }
            list(Obs=sum(TB.L),
                Prob=TB.L/sum(TB.L),
                Pred.class=Pred.class,Dev=Dev.leaf(Lcut[[z]][[1]]),
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
                    kfold=K,oldspV=splvar),
                right=clusterApplyLB(clname,1:length(Rcut),
                    splitnode_par,Rcut,Rtest,Sval=Sval,minSplit=minSplit,
                    dmin=dmin,minBucket=minBucket,top=topN,meth=meth,
                    topnmeth=topn.meth,levelN=levelN,ll=lev,stt=st,Tol=tol,
                    kfold=K,oldspV=splvar)
                )
        }
        else {
            ylev <- table(m[1])
            Pred.class <- names(sort(ylev,decreasing=TRUE)[1])
            Dev <- Dev.leaf(m[[1]])
            fit.tr <- sum(Pred.class == m[[1]])
            if(length(test.d[[1]]) > 0) {
                Dev.test <- Dev.leaf(test.d[[1]])
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
        if (length(Lcut) > 0 && length(Rcut) > 0 ) {
            m<-test.d<-0
            E <-list(split=Sval,
            left=lapply(1:length(Lcut),
                    function(z) {
                    if (length(Lcut[[z]][[1]]) > minSplit && levelN > lev)
                    {
                    sp(Lcut[[z]][[1]],
                        Lcut[[z]],
                        if(length(Ltest) > 0) Ltest[[z]],
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
                        K=K,oldspvar=splvar[z])
                    } else {
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
        right=lapply(1:length(Rcut),
            function(k) {
            if (length(Rcut[[k]][[1]]) > minSplit && levelN > lev )
            {
                sp(Rcut[[k]][[1]],
                Rcut[[k]],
                if(length(Rtest) > 0) Rtest[[k]],
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
                        K=K,oldspvar=splvar[k])
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
    } else {
            ylev <- .Call("tw_table",as.integer(rsp),levels(rsp),PACKAGE="TWIX")
            Pred.class <-attr(which.max(ylev),"names")
            Dev <- Dev.leaf(rsp)
            fit.tr <- sum(Pred.class == rsp)
            if(length(test.d[[1]]) > 0) {
                Dev.test <- Dev.leaf(test.d[[1]])
                id.test <- Pred.class == test.d[[1]]
                fit.test <- .Call("tw_table",as.integer(test.d[[1]][id.test]),
                                levels(test.d[[1]][id.test]),PACKAGE="TWIX")
            } else {
                Dev.test <- fit.test <- 0
                }
            m<-test.d<-0
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
    b <- function(n,Svar,Sp,d,d.t,ftr,ftest,node.cl,sdtr,Obsn,dis,kst,L,R){
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
                    tree[[t<-t+1]]<-list(id=id,dev=dev,dev.test=dev.test,
                        fit.tr=fit.tr,fit.test=fit.test,Splitvar=Splitvar,
                        Splitp=Splitp,Pred.class=node.cl,sd.tr=sd.tr,Obs=Obs,dist=dist,ks.t=ks.t)
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
                        Splitp=Splitp,Pred.class=node.cl,sd.tr=sd.tr,Obs=Obs,dist=dist,ks.t=ks.t)
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
                    tree[[t<-t+1]]<-list(id=id,dev=dev,dev.test=dev.test,
                        fit.tr=fit.tr,fit.test=fit.test,Splitvar=Splitvar,
                        Splitp=Splitp,Pred.class=node.cl,sd.tr=sd.tr,Obs=Obs,dist=dist,ks.t=ks.t)
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
                    if(length(a$left[[k]]) == 3)
                        s(a$left[[k]])
                    else
                        list(list(id=0,Splitvar=0,Splitp=0,dev=0,
                        dev.test=a$left[[k]]$Dev.test,
                        fit.tr=a$left[[k]]$fit.tr,
                        fit.test=a$left[[k]]$fit.test,
                        Pred.class=a$left[[k]]$Pred.class,
                        sd.tr=max(a$left[[k]]$Prob),Obs=a$left[[k]]$Obs,dist=0,ks.t=0)),
                    if(length(a$right[[k]]) == 3)
                        s(a$right[[k]])
                    else
                        list(list(id=0,Splitvar=0,Splitp=0,dev=0,
                        dev.test=a$right[[k]]$Dev.test,
                        fit.tr=a$right[[k]]$fit.tr,
                        fit.test=a$right[[k]]$fit.test,
                        Pred.class=a$right[[k]]$Pred.class,
                        sd.tr=max(a$right[[k]]$Prob),Obs=a$right[[k]]$Obs,dist=0,ks.t=0))
                )
            )
        }
    }
    DF<-id.agr<-agr.id<-1
    if(!is.null(K$split)){
        KK <- s(K)
        id.leng <- sapply(KK,function(x)length(x$id))
        id.out<-which(id.leng == 1)
        no.score<-FALSE
        if(length(id.out) > 0 && length(id.out) < length(KK)){
            KK<-KK[-id.out]
            id.leng <- id.leng[-id.out]
        }
        if(length(id.out) == length(id.leng)){
            KK<-list(KK[[1]])
            gr <- list(KK[[1]])
            no.score<-TRUE
            gr.id<-tw.tic<-gr.tic<-1
        }
        v.fit <- function(x,y){
            y-x$fit.test
        }
        tr.fit <- function(x,z){
            z-x$fit.tr
        }
        for(i in 1:length(KK)) {
            KK[[i]]$dev.test<-K$split[[1]]$Dev.test-KK[[i]]$dev.test
            KK[[i]]$fit.tr<-sum(KK[[i]]$fit.tr)/ntr
            KK[[i]]$fit.test<-sum(KK[[i]]$fit.test)/ntest
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
        no.score<-TRUE
        gr.id<-tw.tic<-gr.tic<-agr.id<-1
    }
    if(!no.score){
    gr <- list(KK[[1]])
    gr.l<-length(gr[[1]]$id[gr[[1]]$id==0])
    if(method == "grid") {
        gr <- list(s(sp(rsp,m,test.data,topn=1,meth="deviance"))[[1]])
        gr[[1]]$fit.tr <- gr[[1]]$fit.tr/ntr
        gr[[1]]$fit.test <- gr[[1]]$fit.test/ntest
    }
    cht<-function(x,n=1,levell=3){
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
    max.fit<-gr.tic<-tw.tic<-DF<-0
    form<-function(obj){
        Var<-splitp<-fit<-colour<-vector()
        for(i in 1:length(obj)){
            Var <- c(Var,as.character(obj[[i]][,1]))
            splitp <- c(splitp,obj[[i]][,2])
            colour <- c(colour,obj[[i]][,3])
            fit <- c(fit,obj[[i]][,4])
        }
        fit <- unclass(as.factor(fit))
        fit<-max(fit)-(fit-1)
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
        fit.tr <- sapply(KK,function(x) x$fit.tr)
        dev.tr <- sapply(KK,function(x) x$dev)
        sd.tr <- lapply(KK,function(x) x$sd.tr[x$Obs != 0])
        Nbuck <- sapply(sd.tr,function(x) length(x))
        sd.tr1 <- sapply(sd.tr,function(x) median(x))
        sd.tr <- sapply(sd.tr,function(x) IQR(x))
        sd.tr <-0.6*sd.tr1/(1+max(sd.tr1))+0.3*(1-sd.tr/(1+max(sd.tr)))
        sd.tr<-sd.tr/max(sd.tr)
        Obs <- sapply(KK,function(x) min(x$Obs[x$Obs != 0]))
        d <- c(which(fit.tr < quantile(fit.tr)[2]),which(dev.tr < quantile(dev.tr)[2]))
        a <- boxplot(Nbuck,plot=FALSE)
        if(minbucket >= 8){
            if(length(KK)> 20)
                id.out <- c(which(Nbuck > a$stats[5,]),which(Nbuck <= a$stats[1,]))
            else
                id.out<-rep(FALSE,length(KK))
            id.agr<-0.6*sd.tr*(1-Nbuck/(1+max(Nbuck)))+0.2*fit.tr/max(fit.tr)+0.2*(1-Obs/max(Obs))
            id.agr[id.out]<-0.1*(id.agr[id.out])
            agr.id <-sort(id.agr,index.return=TRUE,decreasing=TRUE)$ix
            Wbuk<-sd.tr*(1-Nbuck/(1+max(Nbuck)))
            Wbuk<-Wbuk/max(Wbuk)
            if(length(KK)> 20)
                Wbuk[id.out]<-(Wbuk[id.out])*0.25
            DF<-0.15*Wbuk+0.6*fit.tr/max(fit.tr)+0.25*(dev.tr/max(dev.tr))
        }
        if(minbucket < 8){
            if(length(KK)> 20)
                id.out <- c(which(Nbuck <= a$stats[1,]))
            else
                id.out<-rep(FALSE,length(KK))
            id.agr<-0.6*sd.tr*(1-Nbuck/(1+max(Nbuck)))+0.2*fit.tr/max(fit.tr)+0.2*(1-Obs/max(Obs))
            id.agr[id.out]<-0.1*(id.agr[id.out])
            agr.id <-sort(id.agr,index.return=TRUE,decreasing=TRUE)$ix
            Wbuk<-sd.tr*(1-Nbuck/(1+max(Nbuck)))
            Wbuk<-Wbuk/max(Wbuk)
            if(length(KK)> 20)
                Wbuk[id.out]<-(Wbuk[id.out])*0.25
            DF<-0.1*Wbuk+0.7*fit.tr/max(fit.tr)+0.2*(dev.tr/max(dev.tr))
            }
        gr.tic <- DF[1] 
        if(length(KK)> 20)
            DF[d]<-0.1*(DF[d])
        klm<-sort(DF,index.return=TRUE,decreasing=TRUE)$ix
        tw.tic <- DF[klm[1]]
        KK <- KK[klm]
        fit.tr<-dev.tr<-klm<-a<-b<-d<-0
    }
    if(score == 2){
        gr.tic <- KK[[1]]$fit.test
        if(is.data.frame(test.data)){
            KK <- KK[sort(sapply(KK,function(x)x$dev),
                index.return=TRUE,decreasing=TRUE)$ix]
            KK<-KK[sort(sapply(KK,function(x)x$fit.test),
                index.return=TRUE,decreasing=TRUE)$ix]
        }
        else {
            KK <- KK[sort(sapply(KK,function(x)x$dev),
                index.return=TRUE,decreasing=TRUE)$ix]
            KK <- KK[sort(sapply(KK,function(x)x$dev.test),
                index.return=TRUE,decreasing=TRUE)$ix]
            KK<-KK[sort(sapply(KK,function(x)x$fit.test),
                index.return=TRUE,decreasing=TRUE)$ix]
            }
        id.agr<-1:length(KK)
        tw.tic <- KK[[1]]$fit.test
    }
    gr.id <- which(sapply(KK,function(x)x$dev) == gr[[1]]$dev)[1]
    }
    if(!bag) {
        cat("n =",length(KK),"\n")
        cat("Deviance and TIC of the best TWIX-tree:","    ",c(KK[[1]]$dev,tw.tic),"\n")
        if(!is.na(gr.id))
            cat(paste("Deviance and TIC of the greedy tree(Nr.",gr.id,"):"," ",sep =""),
                c(gr[[1]]$dev,gr.tic),"\n")
        else
            cat(paste("Deviance and TIC of the greedy tree:  "," ",sep = ""),
                c(gr[[1]]$dev,gr.tic),"\n")
    }
    greedy.tree <- list(greedy.tree=gr,id=if(!is.na(gr.id)) gr.id else 0)
    class(KK) <- class(greedy.tree) <- "id.tree"
    database <- list(formula=formula(terms(formula,data=m)),
        call=call,greedy.tree=greedy.tree,
        multitree=K,trees=KK,Bad.id=icv,agg.id=agr.id,dist=dist,score=sort(DF,decreasing=TRUE),score2=id.agr)
    class(database) <- "TWIX"
    database
}
