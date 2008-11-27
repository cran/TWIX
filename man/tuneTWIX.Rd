\name{tuneTWIX}
\alias{tuneTWIX}
\title{Parameter tuning.}
\description{
This function tunes hyperparameter minbuck.
}
\usage{
tuneTWIX(formula, data = NULL, minbuck = seq(5,20,by=5),
	maxdepth=30, xval = 10, runs = 10, trace.plot = TRUE)
}
\arguments{
\item{formula}{formula of the form \code{y ~ x1 + x2 + \ldots},
    where \code{y} must be a factor and \code{x1,x2,\ldots} are numeric or factor.}
\item{data}{an optional data frame containing the variables in the model.}
\item{minbuck}{the sampling space for parameter minbuck.}
\item{maxdepth}{set the maximum depth of the final tree.}
\item{xval}{number of cross-validations.}
\item{runs}{number of runs.}
\item{trace.plot}{Should trace plot be ploted?}
}
\seealso{
  \code{\link{plot.TWIX}},\code{\link{TWIX}}
}
\examples{
    library(mlbench)
    data(PimaIndiansDiabetes2)
    Pima <- na.omit(PimaIndiansDiabetes2)

    tuneTWIX(diabetes~.,data=Pima,minbuck=c(5,10,15,20,25),runs=5)
}
\keyword{tree}
