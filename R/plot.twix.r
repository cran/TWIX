plot.bootTWIX <- function(x, ...) plot.TWIX(x, ...)

plot.TWIX <- function(x,sq=1:length(x$trees),type="deviance",
            i.plot=FALSE,size=3,freq=TRUE,breaks = "Sturges",pch=par("pch"),...)
{
    pl.type <- charmatch(type,c("deviance", "ccr", "d&c"))
    if (is.na(pl.type))
        stop("type must be one of 'deviance', 'ccr', 'd&c'")
    n <- length(sq)
    v <- v.test <- vector(,length=n)
    if (!i.plot) {
        plotsect <- function(xx,y,a,b,typ=TRUE) {
            if (x$trees[[1]][[xx]] != 0) {
                for (i in sq) {
                    v[i] <<- x$trees[[i]][[y]]
                    v.test[i] <- x$trees[[i]][[xx]]
                }
            if(typ)
                xxlim <-range(round(min(v)*0.85),round(x$multitree[[1]][[1]]$TD*1.05))
            else
                xxlim <- NULL
            plot(v,v.test,pch=pch,xlab=a,ylab=b,
                xlim=xxlim)
            if(inherits(x, "TWIX"))
                points(x$greedy.tree[[1]][[1]][[y]],
                    x$greedy.tree[[1]][[1]][[xx]],col=2,pch=16)
            }
            else {
                    for (i in sq) {
                    v[i] <<- x$trees[[i]][[y]]
                }
            plot(v,type="o",col=4,pch=pch,ylab=b,)
            if(inherits(x, "TWIX"))
                points(x$greedy.tree[[1]][[1]][[y]],
                which(v ==x$greedy.tree[[1]][[1]][[y]])[1],col=2,pch=16)
                }
        }
        if(pl.type == 0 || pl.type == 1) {
            if(x$trees[[1]]$dev.test == 0) {
                warning("test data not supplied")
                v <- sapply(x$trees,function(y) y$dev)
                hist(v,freq=freq,breaks = breaks,
                    main=paste("Deviances of trees","(n=",
                    length(sq),")"),xlab="Training Deviance")
                abline(v=x$greedy.tree[[1]][[1]]$dev,col=2)
            }
            else {
                par(mfrow=c(2,1))
                plotsect(3,2,"Training Deviance","Test Deviance")
                hist(v,freq=freq,breaks=breaks,main=paste("Deviances of trees",
                    "(n=",length(sq),")"),
                    xlab="Training Deviance",
                    xlim=range(round(min(v)*0.85),
                    round(x$multitree[[1]][[1]]$TD*1.05)))
                abline(v=x$greedy.tree[[1]][[1]]$dev,col=2)
                abline(v=x$multitree[[1]][[1]]$TD,col=3)
            }
        }
        else if(pl.type == 3) {
            if(x$trees[[1]]$fit.test == 0)
                stop("test data not supplied")
            par(mfrow=c(1,1))
            xx<-sapply(x$trees,function(y)y$dev.test)
            yy<-sapply(x$trees,function(y)y$fit.test)
            D<-data.frame(table(xx,yy))
            D<-D[D$Freq != 0,]
            Dx<-as.numeric(levels(D$xx))
            Dy<-as.numeric(levels(D$yy))
            v<-vv<-vector()
            for(i in 1:length(D$x)){
                v[i]<-Dx[D$x[i]]
                vv[i]<-Dy[D$y[i]]
            }
            cex.size<-size * sqrt(D$Freq/max(D$Freq))
            a.size<-range(c(v,vv))
            plot(v,vv,cex=cex.size,pch=16,xlab="Test Deviance",
                ylab="CCR of test data",
                main=paste("Test Deviance vs. Test CCR of",length(xx),"trees"))
            a<-which(round(vv,6) == round(x$greedy.tree[[1]][[1]]$fit.test,6))
            b<-which(round(v,6) == round(x$greedy.tree[[1]][[1]]$dev.test,6))
            id<-a[which(1 == sapply(a,function(x) sum(x==b)))]
            points(v[id],vv[id],cex=cex.size[id],pch=16,col=2)
        }
        if(pl.type == 2) {
            if(x$trees[[1]]$fit.test == 0) stop("test data not supplied")
            par(mfrow=c(1,1))
            xx<-sapply(x$trees,function(y)y$fit.tr)
            yy<-sapply(x$trees,function(y)y$fit.test)
            D<-data.frame(table(xx,yy))
            D<-D[D$Freq != 0,]
            Dx<-as.numeric(levels(D$xx))
            Dy<-as.numeric(levels(D$yy))
            v<-vv<-vector()
            for(i in 1:length(D$x)){
                v[i]<-Dx[D$x[i]]
                vv[i]<-Dy[D$y[i]]
            }
            cex.size<-size * sqrt(D$Freq/max(D$Freq))
            a.size<-range(c(v,vv))
            plot(v,vv,cex=cex.size,pch=16,xlab="CCR of training data",
                ylab="CCR of test data",xlim=a.size,ylim=a.size,
                main=paste("CCR of",length(xx),"trees"))
            a<-which(round(vv,6) == round(x$greedy.tree[[1]][[1]]$fit.test,6))
            b<-which(round(v,6) == round(x$greedy.tree[[1]][[1]]$fit.tr,6))
            id<-a[which(1 == sapply(a,function(x) sum(x==b)))]
            points(v[id],vv[id],cex=cex.size[id],pch=16,col=2)
        }
        par(mfrow=c(1,1))
    }
    else {
        if (length(x$trees[[1]]$dev.test) != 0) {
            f <- f.test <- vector(,length=n)
            for (i in sq) {
                v[i] <- x$trees[[i]]$dev
                v.test[i] <- x$trees[[i]]$dev.test
                f[i] <- x$trees[[i]]$fit.tr
                f.test[i] <- x$trees[[i]]$fit.test
            }
            iplot(v,v.test,ylab="Test Deviance",xlab="Training Deviance")
            idg <- which(v == x$greedy.tree[[1]][[1]]$dev)[1]
            iset.select(idg)
            ig <- rep(0,n)
            ig[idg] <- 2
            iplot(f,f.test,ylab="CCR of test data",xlab="CCR of training data")
            iset.brush(ig)
        if(pl.type == 3)
            iplot(v.test,f.test,xlab="Test Deviance",ylab="CCR of test data")
        }
    }
}
