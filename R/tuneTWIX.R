tuneTWIX<-function(formula, data = NULL, minbuck = seq(5,30,by=5), xval = 10, runs = 10, trace.plot = TRUE){
    call <-match.call()
    m <- match.call(expand=FALSE)
    m$xval <- m$runs <- m$data <- m$minbuck <- NULL
    m <- model.frame(formula,data)
    rsp <- model.extract(m, "response")
    n<-nrow(data)
    runs_in<-length(minbuck)
    CV <- max.er <- rep(0,runs_in)
    CVL <-array(0,c(xval,runs_in,runs))
    xgr <- 1:xval
    for(k in 1:runs){
        id <- sample(rep(xgr, length = n), n)
        for(i in 1:runs_in) {
            for(j in xgr) { 
                test <- id == j 
                train <- !test 
                rpt <- rpart(formula, data=m[train,], 
                    control=rpart.control(cp=0.0001,minbucket=minbuck[i],minsplit=3*minbuck[i],xval=0))
                tiefe <- floor(log2(max(as.numeric(rownames(rpt$frame[rpt$frame[,1] != "<leaf>",])))))+1
                CVL[j,i,k] <- tiefe 
                conf <- sum(predict(rpt,newdata=m[test,], type="class") != m[test,1]) 
                CV[i] <- CV[i] + conf
                max.er[i] <- max.er[i] + nrow(m[test,])
            } 
        }
    }
    CV<-CV/max.er
    ll<-lowess(CV ~ minbuck)
    if(trace.plot){
        plot(minbuck,CV,ylab="error",xlab="minbucket", pch=8,col=4) 
        lines(ll,col=2)
    }
    ij<-which.min(ll$y)
    CVL<-apply(CVL,2,median)
    ji<-which(ll$x[ij] == minbuck)
    out<-data.frame(minsplit=3*ll$x[ij],maxdepth=round(CVL[ji]))
    out
}
