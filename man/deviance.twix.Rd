\name{deviance.TWIX}
\alias{deviance.TWIX}
\alias{deviance.bootTWIX}
\title{Deviance of the TWIX-trees}
\description{
Returns the deviance of a fitted model object.
}
\usage{
\method{deviance}{TWIX}(object,type="training", ...)}
\arguments{
\item{object}{ an object of class \emph{TWIX} or \emph{bootTWIX}
                 for which the deviance is desired.}
\item{type}{ a deviance from test, training or from both data.}
\item{...}{additional arguments.}
}
\value{
The value of the deviance extracted from the object \code{object}.
}
\seealso{
  \code{\link{TWIX}},  \code{\link{predict.TWIX}},  \code{\link{bootTWIX}}
}
\examples{
    data(olives)                                                   
    Tree <- TWIX(Region~.,data=olives,topN=c(5,3),method="local")  
    deviance(Tree)
}
\keyword{tree}
