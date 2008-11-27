\name{trace.plot}
\alias{trace.plot}
\title{Trace plot}
\description{
 A Trace plot shows the structure of the trees
 }
\usage{
trace.plot(obj, sq = 1, quality = NULL, 
		color.palette = topo.colors, alpha = 1) 
}
\arguments{
\item{obj}{formula of the form \code{y ~ x1 + x2 + \ldots},
    where \code{y} must be a factor and \code{x1,x2,\ldots} are numeric or factor.}
\item{sq}{an optional data frame containing the variables in the model.}
\item{quality}{the type of plot: barplot or lines.}
\item{color.palette}{the colors for lines or Bar's.}
\item{alpha}{the alpha transparency, a number in \code{[0,1]}.}
}
\references{
Urbanek, S. (2002). KLIMT - A COSADA Software Project.\cr 
\url{http://rosuda.org/KLIMT}

Urbanek, S. (2005). Following Traces of Lost Models.\cr
\emph{Proc. of the joint Statistical Meetings 2005, Section on Statistical Graphics}\cr
Mira Digital Publishing.
}
\seealso{
  \code{\link{plot.TWIX}},\code{\link{TWIX}},\code{\link{scree.plot}}
}
\examples{
    data(olives)
    Tree1 <- TWIX(Region~.,data=olives[,1:9],topN=c(2,2),method="local")
    get.tree(Tree1)

    ### Trace plot of classification tree
    trace.plot(Tree1,sq=8)

    ### Trace plot of eight classification trees
    trace.plot(Tree1,sq=1:8)
}
\keyword{tree}
