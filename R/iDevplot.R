iDevplot<-function (rsp, data,col=1, ...) {
    require(iplots)
    if (!is.factor(rsp))
        warning("response must be a factor!")
    if (!is.data.frame(data))
        data <- as.data.frame(data)
    data <- na.omit(cbind(rsp, data))
    name <- names(data)
    is <- spoint <-list()
    j<-1
    isl<-dl<-NULL
    w<-iwindow()
    gr1<-igraphics(width=600, height=200,)
    Devplot(data[,1],data[,2],col=col)
    add(w,gr1)
    nw<-dev.list()
    if(is.null(nw))
        idwin<-0
    else 
        idwin<-nw[length(nw)]
    JavaGD(name="Left-node",w=500, h=400)
    JavaGD(name="Right-node",w=500, h=400)
    dev.set(nw)
    g<-igroup(window=w)
    ibutton(text = "Close",handler=function(...) {dev.off(idwin+1);dev.off(idwin+2);dispose(w)}, window=g)
    kl<-icheckbox(text = "show Classes",handler=function(h,...){ get.value(h$obj)} , window=g)
    lab1<-ilabel(" ",window=g)
    update(w)
    for(i in 1:(length(data)-1)){
        F<-splitt(data[,i+1], data[,1], test = TRUE)
        spoint[[i]] <- which(F$dev == max(F$dev)) 
        is[[i]]<-islider(min=1, max=length(F$which),value=spoint[[i]],window=w,
            handler=function(h,x=data,m=i,idw=idwin,points=get.value(kl),...){
                if(is.null(dl)){
                    id <- m+1
                }
                else{
                    id<-which(names(x) == get.value(dl))
                }
                xl <- splitt(x[,id], x[,1], test = TRUE)
                plot(xl$which,xl$dev,type="s",ylab="",xlab="",col=col)
                abline(v=xl$which[get.value(h$obj)],col=2)
                Lcut<-subset(x,x[,id] < xl$which[get.value(h$obj)])
                Rcut<-subset(x,x[,id] >= xl$which[get.value(h$obj)])
                dev.set(idw+1)
                Devplot(Lcut[,1],Lcut[,-1],interactiv=TRUE,col=col,classes=points)
                dev.set(idw+2)
                Devplot(Rcut[,1],Rcut[,-1],interactiv=TRUE,col=col,classes=points)
                dev.set(idw)
            })
        visible(is[[i]],FALSE)
    }
    visible(is[[1]],FALSE)
    ilabel(" Variable ",window=w)
    k<-1
    dl<-idroplist(name[-1], TRUE,window=w,
          handler=function(hh,x=data,obj=is,...){
            id<-which(names(x) == get.value(hh$obj))-1
            if(!visible(obj[[id]])){
                visible(obj[[k]],FALSE)
                visible(obj[[id]],TRUE)
                k<<-id
                }
            Devplot(x[,1],x[,id+1],col=col)                
            })
    visible(w,TRUE)
}
