tuneTWIX<-function(formula, data = NULL, minbuck = seq(5,20,by=5), xval = 10, runs = 10, trace.plot = TRUE){
    call <-match.call()
    m <- match.call(expand=FALSE)
    m$xval <- m$runs <- m$data <- m$minsplit <- NULL
    m <- model.frame(formula,data)
    rsp <- model.extract(m, "response")
    n<-nrow(data)
    minbuck<-rev(minbuck)
    runs_in<-length(minbuck)
    CV <- array(0,c(runs,runs_in))#max.er <- rep(0,runs_in)
    CVL <-array(0,c(xval,runs_in,runs))
    CVN <-array(0,c(xval,runs_in,runs))
    xgr <- 1:xval
    for(k in 1:runs){
        id <- sample(rep(xgr, length = n), n)
        for(i in 1:runs_in) {
            conf<-rep(0,length(xgr))
            for(j in xgr) { 
                test <- id == j 
                train <- !test 
                rpt <- rpart(formula, data=m[train,],parms=list(split="information"), 
                    control=rpart.control(cp=0.0001,minbucket=minbuck[i],minsplit=minbuck[i],
                            xval=10))#,maxcompete=0, maxsurrogate=0, usesurrogate=0))
                tiefe <- floor(log2(max(as.numeric(rownames(rpt$frame[rpt$frame[,1] != "<leaf>",])))))+1
                #maxnode <- length(rpt$frame[rpt$frame[,1] != "<leaf>",1])
                maxnode <- table(predict(rpt,newdata=m[train,], type="class"),m[train,1])
                maxnode <- maxnode[2,1]/length(train)
                #Tdev<-rpt$split[rpt$split[,1] != 0,]
                #CVL[j,i,k] <- sum(Tdev[seq(1,nrow(Tdev),5),3])
                CVL[j,i,k] <- tiefe
                CVN[j,i,k] <- maxnode
                conf[j] <- sum(predict(rpt,newdata=m[test,], type="class") != m[test,1])/length(m[test,1]) 
                #CV[k,i] <- CV[k,i] + conf
                #max.er[i] <- max.er[i] + nrow(m[test,])
            }
            CV[k,i] <- mean(conf)
        }
    }
    CV<-apply(CV,2,mean)#/max.er
    ll<-lowess(CV ~ minbuck)
    if(trace.plot){
        plot(minbuck,CV,ylab="error",xlab="minbuck", pch=8,col=4) 
        lines(ll,col=2)
    }
    #ij<-which.min(ll$y)
    CVL<-apply(CVL,2,median)
    CVN<-apply(CVN,2,mean)
    #ji<-which(ll$x[ij] == minbuck)
    out<-data.frame(minbuck=minbuck[which.min(CV)]#,
                    #maxnode=round(CVN[which.min(CV)],3),
                    #maxdepth=round(CVL[which.min(CV)])
                    )
    out
}
