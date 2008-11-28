\name{Devplot}
\alias{Devplot}
\title{Deviance plot}
\description{
 Deviance plot.
}
\usage{
Devplot(rsp, x, col=1, classes=FALSE, 
	pch=16, nsample=0, ...)
}
\arguments{
\item{rsp}{response variable.}
\item{x}{a dataframe of predictor variables.}
\item{col}{the color for points.}
\item{classes}{Scatterplot of Classes.}
\item{pch}{a vector of plotting characters or symbols.}
\item{nsample}{the number of bootstrap samples without replacement.}
\item{...}{other parameters to be passed through to plotting functions.}
}
\seealso{
  \code{\link{plot.TWIX}},  \code{\link{TWIX}},
  \code{\link{trace.plot}},  \code{\link{scree.plot}}
}
\examples{
    data(olives)

    ### for one predictor variable
    Devplot(olives$Region,olives$oleic)

    ### with classes
    Devplot(olives$Region,olives$oleic,classes=TRUE)

    ### for more as one predictor variable
    Devplot(olives$Region,olives[,1:8])
}
\keyword{tree}
