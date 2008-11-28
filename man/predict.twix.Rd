\name{predict.TWIX}
\alias{predict.TWIX}
\alias{predict.bootTWIX}
\title{
  Predictions from a TWIX object
}
\description{
  The result is a data frame, whose rows are prediction values from 
  appointed tree(s).
}
\usage{
\method{predict}{TWIX}(object,newdata,sq=1,ccr=FALSE,type="class", ...)
}
\arguments{
    \item{object}{an object returned from \code{TWIX} function.}
    \item{newdata}{data frame containing the new data(test data).}
    \item{sq}{Integer vector indicating for which trees predictions are required.}
    \item{ccr}{logical. If \code{TRUE} the result is a list of two components: 
        a data frame with prediction values and correct classification rate of trees.}
    \item{type}{character string indicating the type of predicted value returned. 
        Either \code{class} predicted classes or \code{prob} estimated class 
        probabilities are returned.}
    \item{...}{additional arguments.}
}
\seealso{
\code{\link{bagg}},  \code{\link{TWIX}},  \code{\link{plot.TWIX}}
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

    ### classification 
    pred <- predict(M1,newdata=test,sq=1:2)
    pred

    ### for correct classification rate
    predict(M1,newdata=test,sq=1:2,ccr=TRUE)$CCR

    ### estimated class probabilities
    predict(M1,newdata=test,sq=1,type="prob")
}
\keyword{tree}
