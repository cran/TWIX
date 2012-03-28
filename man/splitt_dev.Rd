\name{splitt}
\alias{splitt}
\title{Internal function of TWIX.}
\description{
This function compute split-points and corresponding deviance gain.}
\usage{
splitt(sv, rsp, meth = "deviance", topn = 1, 
        topn.meth = "complete", lstep = 1,
        test = FALSE, K = 0, level = 0, minbuck = 1)
}
\arguments{                                                                           
\item{sv}{a numeric vector of predicted variable.}
\item{rsp}{ response variable.}
\item{meth}{Which split points will be used?\cr
This can be \code{"deviance"} (default), \code{"grid"} or \code{"local"}.}
\item{topn}{a numeric vector. How many splits will be selected and at which
    level?}
\item{topn.meth}{one of \code{"complete"}(default) or \code{"single"}.}
\item{lstep}{step parameter for method \code{"grid"}.}
\item{test}{parameter for \code{Devplot} funcion.}
\item{K}{k-fold cross-validation.}
\item{level}{the maximum depth of the TWIX tree's.}
\item{minbuck}{the minimum number of observations.}
}
\keyword{tree}
