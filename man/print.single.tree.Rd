\name{print.single.tree}
\alias{print.single.tree}
\title{
  Print tree from single.tree Object.
}
\description{
    This is a method for the generic print() function for objects
    generating by  the function \code{get.tree}.}
\usage{
\method{print}{single.tree}(x, klimt=FALSE, Data=NULL, file="FromR.tree", ...)
}
\arguments{
  \item{x}{an object of class \code{single.tree}.}
  \item{klimt}{logical. If \code{TRUE}, Klimt will be started with
  the tree \code{baum} und dataset \code{Data}.}
  \item{Data}{a data frame. It can be test data or training data.
  This parameter is ignored if \code{klimt == "FALSE"}.}
  \item{file}{a character string naming a file.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\seealso{
  \code{\link{get.tree}},\code{\link{TWIX}}
}
\examples{
  data(olives)
  #Tree<-TWIX(Area~.,data=olives,topN=c(2,2),method="local")
  #Tree1<-get.tree(Tree,n=1)
  #Tree1
  #for Klimt
  #print(Tree1,klimt=TRUE,Data=olives)
}
\keyword{tree}
