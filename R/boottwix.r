bootTWIX<- function(formula, data=NULL,test.data=0,N=1,topN=1,subset=NULL,
                    method="deviance",topn.method="complete",
                    cluster=NULL,minsplit=30,minbucket=round(minsplit/3),
                    Devmin=0.1,level=20,score=1,tol=0.15,splitf="entropy")
{
  BG <-vector("list",length=N)
  if (is.null(cluster)) {
    for(j in 1:N) {
        BGdata <- data[sample(1:nrow(data),nrow(data),replace=TRUE),]
        bg <-TWIX(formula=formula, data=BGdata,test.data,topN=topN,
                    subset=NULL,minsplit=minsplit,minbucket=minbucket,
                    Devmin=Devmin,level=level,bag=TRUE,score=score,
                    method=method,topn.method=topn.method,tol=tol,splitf=splitf)
        BG[[j]] <- bg
        }
  }
  K <- list()
  KK <- lapply(BG,function(x) x$trees)
  for(i in 1:N){
    K <-c(K,list(KK[[i]][[1]]))
    }
  KK<-K
  multitree <- lapply(BG,function(y) y$multitree)
  class(KK) <- "id.tree"
  database <-list(formula=BG[[1]]$formula,call=BG[[1]]$call,
                    multitree=multitree,trees=KK)
  class(database) <- c("bootTWIX",splitf)
  database
}
 
