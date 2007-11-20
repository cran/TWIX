splitt <- function(sv,rsp,svrks=fullrks(sv),
                meth="deviance",topn=1,topn.meth="complete",
                lstep=1,test=FALSE,K=0,level=0)
{
    if(K >= 0.5 || K < 0)
        K <- 0
    n<-length(rsp)
    isCat <- if(is.factor(rsp)) FALSE else TRUE
    if (isCat) {
        stop("\n   Response must be a factor!! \n")
    }
    else{
        rsp <- unclass(rsp)
        }
    if(svrks == 0 && is.factor(sv)) {
        if(topn == 0)
            topn<-32
        Dev<-Which<-vector()
        rtot<-n
        ltot<-0
        sv<-factor(sv)
        t.tot<-table(sv,rsp)
        right<-table(rsp)
        numclass<-length(right)
        numcat<-length(table(sv))
        Wn<-(2**(numcat-1))-1
        tsplit<-gray<-rep(1,length=numcat)  # Right -> 1  Left -> 0
        TD<-sum(-rtot*sapply(right/rtot,nlogn))        
        left<-rep(0,length=numclass)
        Cat_split<-.Call("split_cat", as.matrix(t.tot),
                        as.numeric(TD),
                        as.integer(right),
                        as.integer(numcat),
                        as.integer(numclass),
                        as.numeric(n),
                        as.integer(Wn),
                        PACKAGE="TWIX")
        anf<-1
        enb<-numcat
        Perm<-array(NA,c(Wn,numcat))
        if(length(Cat_split[[3]]) > 0)
        for(i in 1:Wn){
            Perm[i,]<-Cat_split[[3]][anf:enb]
            anf<-enb+1
            enb<-enb+numcat
        }
        if(K != 0)
            K<-trunc(n/(n*K))
        else
            K<-1
        if(test || length(Cat_split[[1]]) <= topn){
            list(dev=Cat_split[[1]],globD=TD,which=Perm,score=rep(K,Wn))
        }
        else {
            id<-sort(Cat_split[[1]],index.return=TRUE,decreasing=TRUE)$ix
            list(dev=Cat_split[[1]][id[1:(topn+1)]],
                globD=TD,
                which=Perm[id[1:(topn+1)], ],
                score=rep(K,length(1:(topn+1))))
        }
    }
    else {
        typ <- charmatch(meth,c("deviance", "local","grid"))
        if (typ == 1)
            meth <- 0
        else if (typ == 2)
            meth <- 1
        else if (typ == 3)
            meth <- 2
# Stichprobe genügend gross =>  cross-valid.
        if( K !=0 && level < 2 && n > 60 && meth != 2){
            csplit<-list()
            xval <- trunc(n/(n*K))
            xgr <- 1:xval
            s <- sample(rep(xgr,length=n),n)
            NN <- n-table(s)
            split_end <- .Call("split_cross",
                        as.numeric(sv),
                        as.integer(rsp),
                        as.integer(NN),
                        as.integer(svrks),
                        as.integer(s),
                        as.integer(xval),
                        PACKAGE="TWIX")
        if (meth == 1) {
            if (length(split_end[[1]]) > 3) {
                dd <- ww <- vector()
                m <- length(split_end[[1]])
                id.l <- 1+
                    which((split_end[[1]][2:(m-1)] > split_end[[1]][1:(m-2)]) &
                        (split_end[[1]][2:(m-1)] > split_end[[1]][3:m]))
                if(sum(split_end[[1]][1] > split_end[[1]]) >= m)
                    id.l <- c(id.l,2)
                if(sum(split_end[[1]][m] > split_end[[1]]) >= m)
                    id.l <- c(id.l,(m+1))
                dd <- split_end[[1]][id.l]
                sc <- split_end[[3]][id.l]
                ww <- split_end[[2]][id.l]
                d <- list(x=dd);wh<-ww
            }
            else {
                d <- list(x=split_end[[1]]);wh<-split_end[[2]];sc<-split_end[[3]]
                }
        }
        else if (meth == 0 && lstep > 1) {
            m <- length(split_end[[1]])
            if (m != 0 ) {
                for (i in 1:3) {
                    split_end[[i]]<-split_end[[i]][seq(1,m,lstep)]
                }
            }
            d<-list(x=split_end[[1]]);wh<-split_end[[2]];sc<-split_end[[3]]
        }
        else {
            d<-list(x=split_end[[1]])
            wh<-split_end[[2]]
            sc<-split_end[[3]]
            }
        if(test || length(d$x) <= topn){
            list(dev=d$x,globD=split_end[[4]],which=wh,score=sc)
        }
        else {
            score<-0.7*d$x/max(d$x) + 0.3*sc/max(sc)
            id <- sort(score,decreasing=TRUE,index.return = TRUE)$ix
            if(topn != 0){
                d$x <- na.omit(d$x[id[1:(topn+1)]])
                wh <- na.omit(wh[id[1:(topn+1)]])
                sc <- na.omit(sc[id[1:(topn+1)]])
            }
            else{
                d$x <- na.omit(d$x[id])
                wh <- na.omit(wh[id])
                sc <- na.omit(sc[id])
            }
            list(dev=d$x,globD=split_end[[4]],which=wh,score=sc)
            }
        }
        else if(meth == 2){
            d <- rep(1,n)
            wh <- seq(min(sv),max(sv),length.out=n)
            if(topn > n)
                topn <- n
            id <- seq(0,n,by=n/(topn+1))
            n_id <- length(id)
            if(topn < n) {
                d <- list(x=d[id[2:(n_id-1)]])
                wh <- wh[id[2:(n_id-1)]]
            }
            else {
                d <- list(x=d[id[2:n_id]])
                wh <- wh[id[2:n_id]]
            }
            list(dev=d$x,globD=10,which=wh)
        }
        else {
            split_end <- .Call("split_single",
                        as.numeric(sv),
                        as.integer(rsp),
                        as.integer(n),
                        as.integer(svrks),
                        PACKAGE="TWIX")
            if (meth == 1) {
                m <- length(split_end[[1]])
                if (m > 3) {
                    dd <- ww <- vector()
                    id.l<-1+which((split_end[[1]][2:(m-1)]>split_end[[1]][1:(m-2)]) &
                            (split_end[[1]][2:(m-1)] > split_end[[1]][3:m]))
                    if(sum(split_end[[1]][1] > split_end[[1]])+1 >= m){
                        id.l <- c(1,id.l)
                        }
                    if(sum(split_end[[1]][m] > split_end[[1]])+1 >= m){
                        id.l <- c(m,id.l)
                        }
                    dd <- split_end[[1]][id.l]
                    ww <- split_end[[2]][id.l]
                    d <- list(x=dd);wh<-ww
                }
                else {
                    d <- list(x=split_end[[1]]);wh<-split_end[[2]]
                }
            }
            else if (meth == 0 && lstep > 1) {
                m <- length(split_end[[1]])
                if (m != 0 ) {
                    for (i in 1:2) {
                        split_end[[i]]<-split_end[[i]][seq(1,m,lstep)]
                    }
                }
                d<-list(x=split_end[[1]]);wh<-split_end[[2]]
            }
            else {
                d<-list(x=split_end[[1]]);wh<-split_end[[2]]
            }
            if(topn.meth == "single" && length(d$x) > (topn+1)) {
                ld <- length(d$x)
                id <- sort(d$x,index.return = TRUE)$ix
                id <- id[ld:(ld-(topn-1))]
                d <- list(x=d$x[id])
                wh <- wh[id]
            }
            if(test || length(d$x) <= topn){
                list(dev=d$x,globD=split_end[[3]],which=wh)
            }
            else {
                id<-sort(d$x,index.return=TRUE,decreasing=TRUE)$ix
                if(topn == 0)
                    list(dev=d$x[id],globD=split_end[[3]],which=wh[id])
                else
                    list(dev=d$x[id[1:(topn+1)]],globD=split_end[[3]],which=wh[id[1:(topn+1)]])
            }
        }
    }
}
#clogn <- function(x) {
#   .C("nlog", as.double(x),PACKAGE="toptree")[[1]]
#}
