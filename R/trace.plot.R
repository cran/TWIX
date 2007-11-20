trace.plot <-function(obj,sq=1,quality=NULL,color.palette=topo.colors,alpha = 1){
        typ <- charmatch(quality,c("deviance", "CCR", "test-CCR", "dev-test"))
        color <- TRUE
        if(is.null(quality)){
            meth <- 4
            color <- FALSE
        }
        else{
            if (typ == 1)
                meth <- 2
            else if (typ == 2)
                meth <- 4
            else if (typ == 3)
                meth <- 5
            else if (typ == 4)
                meth <- 3
            }
    get.splitvar <-function(x,sq=1:length(x$trees),parm="Splitvar",qual=meth) {
        get.ID <- function (m.tr,n=1,id=NULL,which="Splitvar",meth=qual) {
            lev <- j <- i <- k <- m <- 1
            ID<-vector()
            if(class(m.tr)[1] == "TWIX" )
                m.tree <- m.tr$multitree
            else if (class(m.tr)[1] == "bootTWIX")
                m.tree <- m.tr$multitree[[n]]
            if( is.null(id))
                tree.id <- m.tr$trees[[n]]$id
            else
                tree.id <- id
                trase<-array(NA,dim=c(length(tree.id),5))
                ausgabe <- function(m.tree,altSPV,sppold,lev){
                    k<-m
                    root <- i;
                    SP<-m.tree$split[tree.id[root]][[1]]$Splitvar
                    spp<-m.tree$split[tree.id[root]][[1]]$Splitp
                    if(length(spp) > 1){
                        spp<-which(spp == 1)
                    }
                    if(lev > 0){
                        trase[j,]<<-c(altSPV,SP,lev,sppold,spp)
                        j<<-j+1
                    }
                    if (tree.id[i<<-i+1]!=0){
                        m<<-2*k
                        ausgabe(m.tree$left[[tree.id[root]]],SP,spp,lev+1)
                    }
                    if (tree.id[i<<-i+1]!=0){
                        m<<-2*k+1
                        ausgabe(m.tree$right[[tree.id[root]]],SP,spp,lev+1)
                    }
                }
                ausgabe(m.tree,"root",0,1)
                na.omit(trase)
        }
        m<-NULL
        for(i in sq){
            out <- get.ID(x,n=i,which=parm)
            quality <- rep(round(x$trees[[i]][[meth]],4),dim(out)[1])
            m<-rbind(m,cbind(out,quality))
            }
        m
    }
    d.max<-obj$max.range
    xvar<-factor(c(all.vars(obj$formula)[-1]))
    trace.tree<-get.splitvar(obj,sq)

    wh.var<-names(table(as.vector(trace.tree[,1:2])))
    wh.var <- wh.var[-grep("root",wh.var)]

    id.in<-sapply(wh.var,function(x,y){
                        if(x == "root") NA 
                        else which(x==y)}
                ,y=xvar)
    d.max<-d.max[id.in]
    xvar<-wh.var

    x1<-factor(trace.tree[,1],levels=as.character(xvar))
    x2<-factor(trace.tree[,2],levels=as.character(xvar))
    y1<-as.numeric(trace.tree[,3])-1
    y2<-y1+1

    x1.1<-x1;levels(x1.1)<-as.character(d.max)
    x1.2<-x2;levels(x1.2)<-as.character(d.max)

    x1.1<-as.numeric(trace.tree[,4])/as.numeric(as.character(x1.1))
    x1.2<-as.numeric(trace.tree[,5])/as.numeric(as.character(x1.2))
    xx<-seq(1,length(xvar)*2-1,by=2)
    tit.main<-paste("Trace-Plot of",length(sq),"trees",sep=" ")
    
    if(color){
        ############################################################
        ## color legend
        ##
        mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
        w<-(3 + mar.orig[2]) * par("csi") * 2.54
        nf<-layout(matrix(c(2, 1), nc = 2), widths = c(1, lcm(w)))
        #layout.show(nf)
        levels<-pretty(range(as.numeric(trace.tree[,6])), 20)
        col<-color.palette(length(levels)-1)
        par(las = 3)
        mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
        mar <- mar.orig
        mar[4] <- mar[2]
        mar[2] <- 1
        par(mar = mar)
        plot.new()
        plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "r",yaxs = "r")
        rect(0, levels[-length(levels)], 1, levels[-1], col = rev(col))
        axis(4)
        ##
        par(mar = c(5, 2, 4, 2) + 0.1)
        plot(-1,0,xlim=c(1,max(xx)+1),ylim=c(-max(y2)-0.2,0),
                axes = FALSE,xlab="Split-variables",ylab="",main=tit.main)
        colors<-factor(trace.tree[,6])
        levels(colors)<-rev(color.palette(length(levels(colors)),alpha=alpha))
        box()
        for(i in 1:max(y2)){
            l<-1
            for(j in xx){
                lines(c(j,j+1),c(-i,-i))
                text(j-0.1,-i-0.1,xvar[l])
                l<-l+1
            }
        }
        P.Data<-data.frame(x1=xx[unclass(x1)]+x1.1,x2=xx[unclass(x2)]+x1.2,
                            y1=-y1,y2=-y2,col=colors)
        P.Data[is.na(P.Data)] <- (max(xx)+2)/2
        id.out <- which(P.Data[,3] == 0)
        apply(P.Data,1,function(x) lines(x[1:2],x[3:4],col=x[5],lwd=1.5))
        apply(P.Data[-id.out,],1,function(x) points(x[1],x[3],pch="|",cex=1.05))
        apply(P.Data[-id.out,],1,function(x) points(x[2],x[4],pch="|",cex=1.05))
        nf<-layout(matrix(c(1, 1), nc = 1))    
        par(mar = c(5, 4, 4, 2) + 0.1)
    }
    else{
        par(mar = c(5, 2, 4, 2) + 0.1)
        plot(-1,0,xlim=c(1,max(xx)+1),ylim=c(-max(y2)-0.2,0),
                axes = FALSE,xlab="Split-variables",ylab="",main=tit.main)
        box()
        for(i in 1:max(y2)){
            l<-1
            for(j in xx){
                lines(c(j,j+1),c(-i,-i))
                text(j-0.1,-i-0.1,xvar[l])
                l<-l+1
            }
        }
        P.Data<-data.frame(x1=xx[unclass(x1)]+x1.1,x2=xx[unclass(x2)]+x1.2,
                            y1=-y1,y2=-y2,col=rgb(0,0,0,alpha))
        P.Data[is.na(P.Data)] <- (max(xx)+2)/2
        id.out <- which(P.Data[,3] == 0)
        apply(P.Data,1,function(x) lines(x[1:2],x[3:4],col=x[5],lwd=1.5))
        apply(P.Data[-id.out,],1,function(x) points(x[1],x[3],pch="|",cex=1.05))
        apply(P.Data[-id.out,],1,function(x) points(x[2],x[4],pch="|",cex=1.05))
        nf<-layout(matrix(c(1, 1), nc = 1))    
        par(mar = c(5, 4, 4, 2) + 0.1)
    }
}
