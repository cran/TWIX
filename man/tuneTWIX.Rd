\name{tuneTWIX}
\alias{tuneTWIX}
\title{Parameter Tuning.}
\description{
This function tunes hyperparameters minbuck and maxdepth.
}
\usage{
tuneTWIX(formula, data = NULL, minbuck = seq(5, 30, by = 5), 
    xval = 10, runs = 10, trace.plot = TRUE) }
\arguments{
\item{formula}{formula of the form \code{y ~ x1 + x2 + \ldots},
    where \code{y} must be a factor and \code{x1,x2,\ldots} are numeric or factor.}
\item{data}{an optional data frame containing the variables in the model.}
\item{minbuck}{the sampling space for parameter minbuck.}
\item{xval}{number of cross-validations.}
\item{runs}{number of runs.}
\item{trace.plot}{Should trace plot be ploted?}
\item{...}{other parameters to be passed through to plotting functions.}
}
\seealso{
  \code{\link{plot.TWIX}},\code{\link{TWIX}}
}
\examples{
data(olives)
tuneTWIX(Region~.,data=olives[,1:9],minbuck=c(1,5,10,15,20,25),runs=2)
}
\keyword{tree}
