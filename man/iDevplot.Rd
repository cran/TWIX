\name{iDevplot}
\alias{iDevplot}
\title{Interactive Deviance plot}
\description{
 Interactive Deviance plot.
}
\usage{ iDevplot(rsp, data,col=1, ...) } 
\arguments{
\item{rsp}{response variable.}
\item{data}{a dataframe of numeric predictor variables. }
\item{col}{the color for points.}
\item{...}{other parameters to be passed through to plotting functions.}
} \seealso{
  \code{\link{Devplot}},\code{\link{plot.TWIX}},\code{\link{TWIX}}
}
\examples{ 
#data(olives)
#iDevplot(olives$Region,olives[,c(1:8)])
}
\keyword{tree}
