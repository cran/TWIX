\name{Devplot}
\alias{Devplot}
\title{Deviance plot}
\description{
 Deviance plot.
}
\usage{
Devplot(rsp, x,interactiv=FALSE,sample=c(FALSE,0),col=1,classes=FALSE,pch=16, ...)
}
\arguments{
\item{rsp}{response variable.}
\item{x}{a dataframe of predictor variables.}
\item{interactiv}{see \code{\link{iDevplot}}.}
\item{sample}{a dataframe of predictor variables.}
\item{col}{the color for points.}
\item{classes}{Scatterplot of Classes.}
\item{pch}{a vector of plotting characters or symbols.}
\item{...}{other parameters to be passed through to plotting functions.}
}
\seealso{
  \code{\link{iDevplot}},\code{\link{plot.TWIX}},\code{\link{TWIX}}}
}
\examples{
data(olives)
#Devplot(olives$Region,olives[,1:8])
}
\keyword{tree}
