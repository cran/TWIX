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
\method{bagg}{default}(object,newdata = NULL, sq=1:length(object$trees),
        aggregation = "weighted", type = "class", \dots)
\method{bagg}{TWIX}(object,\dots)
\method{bagg}{bootTWIX}(object,\dots)
}
\arguments{
\item{object}{object of classes \code{TWIX} or \code{bootTWIX}.}
\item{newdata}{a data frame of new observations.}
\item{sq}{Integer vector indicating for which trees predictions are required.}
\item{aggregation}{character string specifying how to aggregate. There are 
	three ways to aggregate the TWIX trees. One of them is the class majority voting 
	(\code{aggregation="majority"}), another method is the weighted aggregation 
	(\code{aggregation="weighted"}).}
\item{type}{character string indicating the type of predicted value returned. 
	Either \code{class} predicted classes or \code{prob} estimated class 
	probabilities are returned.}
\item{\dots}{additional arguments.}
}

\seealso{                                             
  \code{\link{TWIX}},  \code{\link{predict.TWIX}},  \code{\link{bootTWIX}}
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

    ### TWIX Ensemble
    M1 <- TWIX(chd~.,data=train,topN=c(9,5),topn.method="single")

    ### Bagging with greedy decision trees as base classifier
    M2 <- bootTWIX(chd~.,data=train,nbagg=50)

    ### Bagging with the p-value adjusted classification trees
    ### as base classifier
    M3 <- bootTWIX(chd~.,data=train,nbagg=50,splitf="p-adj",Devmin=0.01)

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
