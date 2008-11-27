Devplot <- function(rsp, x, col=1, classes=FALSE, pch=16, nsample=0, ...) 
{
    name <- names(x)
    outl <- fact <- FALSE
    kl <- vector()
	interactiv <- FALSE
    if (!is.factor(rsp)) 
        stop("response must be a factor!")
    if (!is.data.frame(x) && !inherits(x, "TWIX")){
        x <- as.data.frame(x)
        fact <-  sapply(x,is.factor)
        if(sum(fact) > 0)
            stop("some predictor variables are factors")
        }
    if(inherits(x, "TWIX")){
        outl <- TRUE
        Lx <- x[[2]]
        pcol <- x[[3]]
        x <- as.data.frame(x[[1]])        
        fact <-  sapply(x,is.factor)
        if(sum(fact) >= ncol(x))
            stop("all predictor variables are factors")
        if(sum(fact) > 0)
            warning(paste("variables","(",paste(names(fact[fact]),collapse=", "),")"," are factors!"))
        kl <- which(fact == TRUE)
        idl <- sapply(names(Lx),function(x,y){which(x == y)},y=names(x))
        outf<-vector()
        for(i in 1:length(kl))
            outf<-c(outf,which(kl[i] == idl))
        if(length(outf) > 0){
            idl <- idl[-outf]
            Lx <-Lx[-outf]
            }
        id.il <- fullrks(idl)
        color <- c(2,4,4)
        color2 <- gray((0:pcol)/(pcol+2))
        x <- na.omit(cbind(rsp, x))
        n <- ncol(x[!fact])        
    }
    else{
        x <- na.omit(cbind(rsp, x))
        n <- ncol(x)
        }
    N <- ncol(x)
    newY <- function(x,y){
        lenax<-round(max(x$dev))
        Ytick<-pretty(x$dev)
        ntick<-length(Ytick)
        diftick<-Ytick[ntick]-Ytick[ntick-1]
        newL<-unclass(y)*(diftick/5)+lenax
        }
    par(mfrow = c(j <- ceiling((n - 1)/2), ceiling((n - 1)/j)),
        mar = c(2.5, 4, 3, 1) + 0.1)
    if (n == 2) {
        D <- splitt(x[, 2], x[, 1], test = TRUE)
        if(interactiv){
            newL<-newY(D,x[,1])
            plot(D$which, D$dev,type="s",xlab = "", main = "", ylab = "",
                cex=0.5,pch=16,col=col,ylim=c(0,max(newL)))
            if(classes){
                points(x[,2],newL,col=unclass(x[,1]),pch=pch, ...)
                }
        }
        else if(nsample > 0){
            out <- list()
            for(m in 1:nsample){
                icv <- sample(1:nrow(x),0.632*nrow(x))
                out[[m]] <- splitt(x[icv,2],x[icv,1],test=TRUE)
                }
            yaxlim <- max(sapply(out,function(x) max(x$dev)))
            plot(out[[1]]$which,out[[1]]$dev,"l",col=col,xlab = "Splitpoints",
                    ylab =expression(Delta*i(s,t)),ylim=c(0,yaxlim),xlim=range(x[,2]))
            for(i in 2:nsample)
				lines(out[[i]]$which,out[[i]]$dev,col=col)
        }
        else {
            newL<-newY(D,x[,1])
            plot(D$which, D$dev,xlab = "Splitpoints", main = "",
                ylab =expression(Delta*i(s,t)),type = "s",lwd=1,ylim=c(0,max(newL)), ...)
            if(classes)
                points(x[,2],newL,col=unclass(x[,1]),pch=pch, ...)
            if(outl){
                apply(Lx[[1]],1,function(x,y)abline(v=x[1],col=y[x[2]]),y=color)
                apply(Lx[[1]],1,function(x,y)points(x[1],0,col=y[x[3]],pch=16),y=color2)
            }
        }
    }
    else {
        D <- lapply(2:N, function(z, y) {
            splitt(y[, z], y[, 1], test = TRUE)
        }, y = x)
        if(interactiv){
            if(classes)
                newL<-lapply(D,newY,y=x[,1])
            else
                newL<-lapply(D,function(x) max(x$dev))
            rangeY<-sapply(newL,max)
            inn <- which.max(rangeY)
            nmax <- rangeY[inn]
            size <- range(unlist(sapply(1:(n - 1), function(y) range(D[[y]]$dev))))
            name <- names(x)
            for (i in 1:(n - 1)){
                plot(D[[i]]$which, D[[i]]$dev, main = name[[i+1]], 
                    xlab = "", ylim = c(0,nmax), ylab = "", type = "s",col=col, ...)
                if(classes)
                    points(x[,i+1],newL[[inn]],cex=0.5,pch=pch,col=unclass(x[,1]))
                }
        }
        else if(nsample > 0){
            D <- lapply(2:N, function(z,y,n) {
                    out <- list()
                    for(g in 1:n){
                        icv <- sample(1:nrow(y),0.632*nrow(y))
                        out[[g]] <- splitt(y[icv, z], y[icv, 1], test = TRUE)
                        }
                    out
                    },y=x,n=nsample)
            yaxlim <- max(unlist(lapply(D,function(x) lapply(x,function(x)(max(x$dev))))))
            xaxmax <- lapply(D,function(x) max(sapply(x,function(x)(max(x$which)))))
            xaxmin <- lapply(D,function(x) min(sapply(x,function(x)(min(x$which)))))
            if(length(kl > 0))
                sql <- (1:(N - 1))[-kl]
            else
                sql <- (1:(N - 1))
            for (i in sql){
                plot(D[[i]][[1]]$which,D[[i]][[1]]$dev,"l",col=col,xlab = "Splitpoints", main = name[[i]],
                    ylab =expression(Delta*i(s,t)),ylim=c(0,yaxlim),xlim=c(xaxmin[[i]],xaxmax[[i]]))
                for(h in 2:nsample)
                    lines(D[[i]][[h]]$which,D[[i]][[h]]$dev,col=col)
                }
        }
        else{
            if(classes)
                newL<-lapply(D,newY,y=x[,1])
            else
                newL<-lapply(D,function(x) max(x$dev))
            rangeY<-sapply(newL,max)
            inn <- which.max(rangeY)
            nmax <- rangeY[inn]
            size <- range(unlist(sapply(1:(n - 1), function(y) range(D[[y]]$dev))))
            name <- names(x)
            p<-1
            if(length(kl > 0))
                sql <- (1:(N - 1))[-kl]
            else
                sql <- (1:(N - 1))
            for (i in sql){
                plot(D[[i]]$which, D[[i]]$dev, main = name[[i+1]], 
                    xlab = "", ylim = c(0,nmax), ylab = "", type = "s",col=col, ...)
                if(classes)
                    points(x[,i+1],newL[[inn]],cex=0.5,pch=pch,col=unclass(x[,1]))
                if(outl && i == idl[id.il[p]] && p <= length(id.il)){
                    xx<-Lx[[id.il[p]]]
                    ij<-fullrks(xx[,1])
                    xx <- xx[ij,]
                    anf <- 0
                    for(j in 1:nrow(Lx[[id.il[p]]])){
                        abline(v=xx[j,1],col=color[xx[j,2]])
                        if(j > 1 && xx[j-1,1] == xx[j,1]){
                            anf <- anf+nmax*0.01
                            points(xx[j,1],anf,col=color2[xx[j,3]],cex=1.2,pch=16)
                        }
                        else{
                            points(xx[j,1],0,col=color2[xx[j,3]],cex=1.2,pch=16)
                            anf <- 0
                            }
                    }
                    p<-p+1
                }
            }
        }
    par(mfrow = c(1, 1))
    }
}
