\name{bagg.default}
\alias{bagg}
\alias{bagg.default}
\alias{bagg.bootTWIX}
\alias{bagg.TWIX}
\title{ Predictions from TWIX's or Bagging Trees }
\description{
  Prediction of a new observation based on multiple trees.  
}
\usage{
\method{bagg}{default}(object,newdata=NULL,sq=1:length(object$trees),
		aggregation="weighted",\dots)
\method{bagg}{TWIX}(object,\dots)
\method{bagg}{bootTWIX}(object,\dots)
}
\arguments{
\item{object}{object of classes \code{TWIX} or \code{bootTWIX}.}
\item{newdata}{a data frame of new observations.}
\item{sq}{Integer vector indicating for which trees predictions are required.}
\item{aggregation}{character string specifying how to aggregate. There are 
	two ways to aggregate the TWIX trees. One of them is the class majority voting 
	(\code{aggregation="majority"}) and another method is the weighted aggregation 
	(\code{aggregation="weighted"}).}
\item{\dots}{additional arguments.}
}

\seealso{                                             
  \code{\link{TWIX}}, \code{\link{predict.TWIX}}
  \code{\link{bootTWIX}}
}
  
\examples{
    library(ElemStatLearn)
    data(SAheart)

    ### response variable must be a factor
    SAheart$chd <- factor(SAheart$chd) 

    ### test and train data
    ###
    set.seed(1234)
    icv <- sample(nrow(SAheart),nrow(SAheart)/3)
    itr <- setdiff(1:nrow(SAheart),icv)
    train <- SAheart[itr,]
    test <- SAheart[icv,]
  
    M1 <- TWIX(chd~.,data=train,topN=c(4,3),topn.method="single")
    M2 <- bootTWIX(chd~.,data=train,nbagg=30)
    M3 <- bootTWIX(chd~.,data=train,nbagg=30,splitf="p-adj",Devmin=0.01)

    preda <- bagg(M1,test,sq=1:length(M1$trees),aggregation="majority")
    predb <- bagg(M1,test,sq=1:length(M1$trees),aggregation="weighted")
    pred1 <- predict(M2,test,sq=1:length(M2$trees))
    pred2 <- predict(M3,test,sq=1:length(M3$trees))

    ###
    ### CCR's
  
    sum(preda == test$chd)/nrow(test)
    sum(predb == test$chd)/nrow(test)
    sum(pred1 == test$chd)/nrow(test)
    sum(pred2 == test$chd)/nrow(test)
} 
\keyword{tree}
