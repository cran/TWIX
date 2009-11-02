\name{print.single.tree}
\alias{print.single.tree}
\title{
  Print tree from single.tree object.
}
\description{
    This is a method for the generic print() function for objects
    generating by the function \code{get.tree}.}
\usage{
\method{print}{single.tree}(x, klimt=FALSE, Data=NULL, file="FromR.tree", ...)
}
\arguments{
  \item{x}{an object of class \code{single.tree}.}
  \item{klimt}{logical. If \code{TRUE}, Klimt will be started with
  the tree \code{baum} und dataset \code{Data}.}
  \item{Data}{a data frame. It can be test or training data.\cr
  This parameter is ignored if \code{klimt == "FALSE"}.}
  \item{file}{a character string naming a file.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\references{
Urbanek Simon (2002). KLIMT - A COSADA Software Project.\cr 
\url{http://rosuda.org/KLIMT}
}
\seealso{
  \code{\link{get.tree}},  \code{\link{TWIX}}
}
\examples{
    data(olives)
    Tree <- TWIX(Area~.,data=olives,topN=c(2,2),method="local")
    Tree <- get.tree(Tree,n=1)
    Tree

    ### for further analysis in KLIMT
    print(Tree,klimt=TRUE,Data=olives)
}
\keyword{print} \keyword{tree}
