\name{predict.TWIX}
\alias{predict.TWIX}
\alias{predict.bootTWIX}
\title{
  Predictions from a TWIX Object
}
\description{
  The result is a data frame, whose rows are prediction values from 
  appointed tree(s).
}
\usage{
\method{predict}{TWIX}(object,newdata,sq=1,ccr=FALSE,...)
}
\arguments{
  \item{object}{an object returned from \code{TWIX} function.}
  \item{newdata}{data frame containing the new data(test data).}
  \item{sq}{Integer vector indicating for which trees predictions
  are required.}
  \item{ccr}{logical. If \code{TRUE} the result is a list of two
  components: a data frame with prediction values and correct 
  classification rate of trees. This parameter can be ignored, if
  the function \code{TWIX} has been called call with test data
  (\code{test.data=test}).}
  \item{...}{additional arguments affecting the predictions produced}
}
\seealso{
  \code{\link{TWIX}}, \code{\link{plot.TWIX}}
}
\examples{
data(olives)
set.seed(1234)
i <- sample(572,150)
ic <- setdiff(1:572,i)
training <- olives[ic,]
test <- olives[i,]

Tree<-TWIX(Region~.,training,topN=c(5,2),method="local")
pred<-predict(Tree,newdata=test,sq=1:2)

#
## for Correct Classification Rate
#

predict(Tree,newdata=test,sq=1:2,ccr=TRUE)$CCR
}
\keyword{tree}
