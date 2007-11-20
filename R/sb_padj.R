    badj <- function(n,Svar,Sp,d,d.t,ftr,ftest,node.cl,sdtr,Obsn,dis,kst,L,R){
        tree <- list()
		t<-0
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
                    tree[[t<-t+1]]<-list(id=id,dev=dev,dev.test=dev.test,
                        fit.tr=fit.tr,fit.test=fit.test,Splitvar=Splitvar,
                        Splitp=Splitp,Pred.class=node.cl,sd.tr=sd.tr,Obs=Obs,dist=dist,ks.t=ks.t)
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
                        Splitp=Splitp,Pred.class=node.cl,sd.tr=sd.tr,Obs=Obs,dist=dist,ks.t=ks.t)
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
                    tree[[t<-t+1]]<-list(id=id,dev=dev,dev.test=dev.test,
                        fit.tr=fit.tr,fit.test=fit.test,Splitvar=Splitvar,
                        Splitp=Splitp,Pred.class=node.cl,sd.tr=sd.tr,Obs=Obs,dist=dist,ks.t=ks.t)
                }
            }
        }
    tree
    }
	
    sadj <- function(a) {
        id <- list()
        for (k in 1:length(a$split)) {
            id <- c(id,
                badj(k,a$split[[k]]$Splitvar,a$split[[k]]$Splitp,a$split[[k]]$Dev,
                        a$split[[k]]$Dev.test,0,0,
                        a$split[[k]]$Pred.class,
                        max(a$split[[k]]$Prob),0,a$split[[k]]$dist,
                        a$split[[k]]$ks.t,
                    if(length(a$left[[k]]) == 3)
                        sadj(a$left[[k]])
                    else
                        list(list(id=0,Splitvar=0,Splitp=0,dev=a$left[[k]]$Dev,
                        dev.test=a$left[[k]]$Dev.test,
                        fit.tr=a$left[[k]]$fit.tr,
                        fit.test=a$left[[k]]$fit.test,
                        Pred.class=a$left[[k]]$Pred.class,
                        sd.tr=max(a$left[[k]]$Prob),Obs=a$left[[k]]$Obs,dist=0,ks.t=0)),
                    if(length(a$right[[k]]) == 3)
                        sadj(a$right[[k]])
                    else
                        list(list(id=0,Splitvar=0,Splitp=0,dev=a$right[[k]]$Dev,
                        dev.test=a$right[[k]]$Dev.test,
                        fit.tr=a$right[[k]]$fit.tr,
                        fit.test=a$right[[k]]$fit.test,
                        Pred.class=a$right[[k]]$Pred.class,
                        sd.tr=max(a$right[[k]]$Prob),Obs=a$right[[k]]$Obs,dist=0,ks.t=0))
                )
            )
        }
    }
