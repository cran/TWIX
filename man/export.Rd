\name{export}
\alias{export}
\title{Export TWIX-trees for KLIMT}
\description{
 Export TWIX-trees for KLIMT.
}
\usage{
export(x, sq = 1, directory = "ForKlimt")
}
\arguments{
\item{x}{an object of class \code{TWIX}.}
\item{sq}{a vector of tree IDs.}
\item{directory}{a name of directory, which will be created.}
}
\seealso{
  \code{\link{predict.TWIX}},\code{\link{print.TWIX}}}
}
\examples{
data(olives)
#Tree<-TWIX(Region~.,data=olives,topN=c(3,2),method="local",level=5)
#export(Tree,c(1,3))
}
\keyword{tree}
