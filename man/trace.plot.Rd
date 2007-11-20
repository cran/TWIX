\name{trace.plot}
\alias{trace.plot}
\title{Trace-plot}
\description{
 A Trace plot shows
 }
\usage{
trace.plot(obj, sq = 1, quality = NULL, color.palette = topo.colors, 
		alpha = 1) 
}
\arguments{
\item{obj}{formula of the form \code{y ~ x1 + x2 + \ldots},
    where \code{y} must be a factor and \code{x1,x2,\ldots} are numeric or factor.}
\item{sq}{an optional data frame containing the variables in the model.}
\item{quality}{the type of plot: barplot or lines.}
\item{color.palette}{The colors for lines or Bar's.}
\item{alpha}{the color for points.}
}
\seealso{
  \code{\link{iDevplot}},\code{\link{plot.TWIX}},\code{\link{TWIX}},\code{\link{scree.plot}}
}
\examples{
data(olives)
Tree1<-TWIX(Region~.,data=olives[,1:9],topN=c(2,2),method="local")
trace.plot(Tree1,sq=8)
}
\keyword{tree}
