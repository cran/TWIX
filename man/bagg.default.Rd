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
\method{bagg}{default}(object,data=NULL,sq=1:10,\dots)
\method{bagg}{TWIX}(object,\dots)
\method{bagg}{bootTWIX}(object,\dots)
}
\arguments{
  \item{object}{object of classes \code{TWIX} or
                \code{bootTWIX}.}
  \item{data}{a data frame of new observations. }
  \item{sq}{Integer vector indicating for which trees predictions are required.}
  \item{\dots}{\dots}
}
\seealso{                                             
  \code{\link{TWIX}}, \code{\link{predict.TWIX}}
  \code{\link{bootTWIX}}
}  
\examples{
data(olives)
#Tree <- TWIX(Region~.,data=olives,topN=c(5,3),method="local")
#Tree1 <- bootTWIX(Region~.,data=olives,topN=c(3,1),N=10)
#pred <- bagg(Tree,olives,sq=1:10)
#pred1 <- bagg(Tree1,olives,sq=1:10)
#
#CCR's
#sum(pred==olives$Region)/nrow(olives)
#sum(pred1==olives$Region)/nrow(olives)
} 
\keyword{tree}
