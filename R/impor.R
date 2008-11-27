impor<-function(newdata,y,runs){
    xval <- 4
	if(is.data.frame(newdata))
		newdata <- data.frame(newdata)
	n <- nrow(newdata)
    CVL <-array(0,dim=c(xval*runs,ncol(newdata)))
    CV <-array(0,dim=c(xval*runs,ncol(newdata)))
    xgr <- 1:xval
    v<-1
    for(k in 1:runs){
        id <- sample(rep(xgr, length.out = n), n)
        for(j in xgr) { 
            test <- id == j 
            train <- !test 
            S<-lapply(newdata[train,],splitt,y[train])
			CV[v,]<-.Call("Dev_oob",S,newdata[test,],as.numeric(y[test]),PACKAGE="TWIX")
            CVL[v,]<-unlist(lapply(S,function(x)(x$dev[1])))
            v<-v+1
        }
    }
    b<-apply(CVL,2,median)
    a<-apply(CV,2,median)
    #cat("Wicht. Var nach train: ",names(data.train)[rev(fullrks(b))[1:3]],"\n")
    #cat("Wicht. Var nach test: ",names(newdata)[rev(fullrks(a))[1:3]],"************* \n")
    #cat("Wicht. Var nach Ts_Tr: ",names(data.train)[rev(fullrks(apply(rbind(a,b),2,mean)))[1:3]],"************* \n")
    #cat("\n ------------------ \n")
    #names(newdata)[rev(fullrks(a))[1:3]]
	list(dev.tr=b,dev.test=a)
}
