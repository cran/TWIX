TWIX <- function(formula, data=NULL, test.data=NULL, subset=NULL,
        method="deviance", topn.method="complete", minsplit=20,
		minbucket=round(minsplit/3), topN=1, splitf="deviance",
		Devmin=0.01, tol=0.25, cp=0.01, level=30, st=1, score=1, k=0,
		cluster=NULL, seed.cluster=NULL, cl.level=1, multicore=FALSE,
		trace=TRUE, trace.plot=FALSE, ...)
{
    call <- match.call()
    m <- match.call(expand.dots = FALSE)
    m$method <- m$topn.method <- m$test.data <- NULL
    m$minsplit <- m$minbucket <- m$trace <- m$cp <- NULL
    m$Devmin <- m$topN <- m$level <- m$st <- m$tol <- m$... <- NULL
    m$score <- m$k <- m$trace.plot <- m$robust<- m$splitf <- NULL
	m$cluster <- m$cl.level <- m$multicore <- NULL
    m[[1]] <- as.name("model.frame.default")
    m <- eval(m, parent.frame())
	
	data.Na <- apply(m, 1, function(x) any(is.na(x)))
	if(sum(data.Na) > 0){
		stop("\n Missing values in training data! \n\n")
		}
    if(splitf == "p-adj")
        Devmin <- 1-Devmin
    if(is.null(test.data) || is.data.frame(test.data)){
        m <- twix.data(m,splitf)
    }
    else if(is.matrix(test.data)){
        stop("\n   test.data must be a data.frame!\n")
        }
    else if(length(test.data) == 1){
        test.data <- NULL
        }
    if(!is.null(test.data)) {
        test.data <- model.frame(formula,test.data)
		test.data <- twix.data(test.data,splitf)
        ntest <- dim(test.data)[1]
    }
	else{
        ntest <- 1
	}
    ntr <- dim(m)[1]
    rsp <- m[[1]]
	robust <- FALSE
    j <- i <- 1
    if(!is.null(names(list(...)))){
		if(names(list(...)) == "robust")
			robust <- TRUE
	}
    if(method=="grid" && topn.method != "single")
        topn.method <- "single"
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
	if(splitf != "p-adj" & splitf != "deviance")
			stop("\n   'splitf' must be 'deviance' or 'p-adj'! \n")
    if(!is.null(cluster)){
		if(!is.null(seed.cluster)){
			clusterSetRNGStream(cluster, iseed=seed.cluster)
		}else{
			clusterSetRNGStream(cluster, iseed=sample(1:9999,6))
		}
        clusterEvalQ(cluster, library(TWIX))
	}
    sp <- function(rsp, m, test.d=test.data, dmin=Devmin, minSplit=minsplit,
                minBucket=minbucket, clname=cluster, cllevel=cl.level, topn=topN,
				topn.meth=topn.method, levelN=level,lev=0, meth=method, lstep=st,
				tl=tol, K=k, rob=robust, splf=splitf,multic=multicore)
    {
    n <- length(m)
    E <- vector("list",3)
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
    if(lev == 2 && meth == 2) {
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
		if(rob && lev < 3){
			imp <- importance(m[,2:n],rsp,runs=10)
			S <- .Call("var_split_dev", m, as.integer(meth), 
						as.integer(lstep), as.integer(n_cut),
						as.logical(FALSE), as.numeric(K),
						as.integer(minBucket), as.integer(lev),
						as.logical(TRUE), as.numeric(tl),
						PACKAGE="TWIX")
			for(v in 1:length(S)){
				a <- ((S[[v]][[1]]/max(S[[v]][[1]])))*imp[[2]][v]*4
				a[is.na(a)] <- 0
				S[[v]][[1]] <- a
			}
			S <- .Call("split_summary_dev",S,as.numeric(tl),PACKAGE="TWIX")
		}
		else{
			S <- .Call("var_split_dev", m, as.integer(meth), 
						as.integer(lstep), as.integer(n_cut),
						as.logical(FALSE), as.numeric(K),
						as.integer(minBucket), as.integer(lev),
						as.logical(FALSE), as.numeric(tl),
						PACKAGE="TWIX")
		}
	}
	else{
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

    splitnode_par <- function(z,Lcut,Ltest,Sval,dmiN,minSplit,minBucket,top,
                            meth,topnmeth,levelN,ll,stt,Tol,kfold,splitfun,Rob)
    {
        if(length(Lcut[[z]][[1]]) >= minSplit && levelN > ll )
        {
            sp.slave(Lcut[[z]][[1]],Lcut[[z]],
                if(length(Ltest) >= z && !is.null(Ltest[[z]])) Ltest[[z]],
                Dmin=dmiN,minsplit=minSplit,
                minbucket=minBucket,topN=top,
                method=meth,topn.method=topnmeth,
                level=levelN,lev=ll,st=stt,tol=Tol,
                K=kfold,splitf=splitfun,robust=Rob)
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
    if(lev >= cllevel && (!is.null(clname) || multic)) {
        if (h_level > 1) {
			if(!multic){
				E <-list(split=Sval,
						 left=clusterApplyLB(clname,1:(h_level-1),
											 splitnode_par,Lcut,Ltest,Sval=Sval,dmiN=dmin,minSplit=minSplit,
											 minBucket=minBucket,top=topN,meth=meth,
											 topnmeth=topn.meth,levelN=levelN,ll=lev,stt=st,Tol=tol,
											 kfold=K,splitfun=splf,rob),
						 right=clusterApplyLB(clname,1:(h_level-1),
											  splitnode_par,Rcut,Rtest,Sval=Sval,dmiN=dmin,minSplit=minSplit,
											  minBucket=minBucket,top=topN,meth=meth,
											  topnmeth=topn.meth,levelN=levelN,ll=lev,stt=st,Tol=tol,
											  kfold=K,splitfun=splf,rob)
						 )
			}else{
				E <-list(split=Sval,
						 left=mclapply(1:(h_level-1),
									splitnode_par,Lcut,Ltest,Sval=Sval,dmiN=dmin,minSplit=minSplit,
									minBucket=minBucket,top=topN,meth=meth,
									topnmeth=topn.meth,levelN=levelN,ll=lev,stt=st,Tol=tol,
									kfold=K,splitfun=splf,Rob=rob, mc.preschedule = FALSE),
						 right=mclapply(1:(h_level-1),
									splitnode_par,Rcut,Rtest,Sval=Sval,dmiN=dmin,minSplit=minSplit,
									minBucket=minBucket,top=topN,meth=meth,
									topnmeth=topn.meth,levelN=levelN,ll=lev,stt=st,Tol=tol,
									kfold=K,splitfun=splf,Rob=rob, mc.preschedule = FALSE)
						 )
			}
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
        if(h_level > 1){
            E <- list(split=Sval,
            left=my.lapply(1:(h_level-1), 
					function(z){
                    if(length(Lcut[[z]][[1]]) >= minSplit && levelN > lev){
                    sp(Lcut[[z]][[1]],
                        Lcut[[z]],
                        if(length(Ltest) >= z && length(Ltest[[z]][[1]]) > 0) Ltest[[z]],
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
                        K=K,rob=rob)
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
			right=my.lapply(1:(h_level-1),
				function(k) {
				if(length(Rcut[[k]][[1]]) >= minSplit && levelN > lev )
				{
                sp(Rcut[[k]][[1]],
					Rcut[[k]],
					if(length(Rtest) >= k && length(Rtest[[k]][[1]]) > 0) Rtest[[k]],
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
					K=K,rob=rob)
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
			#	rm(list=c("test.d","Ltest","Rtest"))
			}
			else{
                Dev.test <- fit.test <- 0
			}
            E <- list(Obs=sum(ylev),Prob=ylev/sum(ylev),
                        Pred.class=Pred.class,Dev=Dev,
                        Dev.test=Dev.test,fit.tr=fit.tr,
                        fit.test=fit.test)
		}
	}
    E
  }
    K <- sp(rsp,m,test.data)
	DF <- id.agr <- agr.id <- 1
    if(!is.null(K$split)){
		if(splitf != "p-adj"){
			KK <- s(K,K$split[[1]]$Obs-K$split[[1]]$fit.tr,cp)
		}
		else{
			KK <- s_padj(K)
		}
        id.leng <- sapply(KK,function(x)length(x$id))
        id.out <- id.leng == 1
        no.score <- FALSE
        if(sum(id.out) > 0 && sum(id.out) < length(KK)){
            KK <- KK[!id.out]
            id.leng <- id.leng[!id.out]
        }
        if(sum(id.out) == length(id.leng)){
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
        K <- list(split=list(K))
		
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
        klm <- my.sort(fit.tr_certaintly,index.return=TRUE,decreasing=TRUE)
        tw.tic <- klm$x[1]
        KK <- KK[klm$ix]
		DF <- klm$x
        fit.tr_certaintly <- klm <- 0
    }
    if(score == 3){
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
    if(score == 2){
        if(is.null(test.data)){
			gr.tic <- KK[[1]]$fit.tr
            KK <- KK[my.sort(sapply(KK,function(x)x$dev),
                index.return=TRUE,decreasing=TRUE)$ix]
            KK<-KK[my.sort(sapply(KK,function(x)x$fit.tr),
                index.return=TRUE,decreasing=TRUE)$ix]
        }
        else {
			gr.tic <- KK[[1]]$fit.test
            KK <- KK[my.sort(sapply(KK,function(x)x$dev),
                index.return=TRUE,decreasing=TRUE)$ix]
            KK <- KK[my.sort(sapply(KK,function(x)x$dev.test),
                index.return=TRUE,decreasing=TRUE)$ix]
            KK<-KK[my.sort(sapply(KK,function(x)x$fit.test),
                index.return=TRUE,decreasing=TRUE)$ix]
            }
		DF <- 1:length(KK)
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
    if(trace) {
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
