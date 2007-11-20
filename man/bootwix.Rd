\name{bootTWIX}
\alias{bootTWIX}
\title{
Bootstrap of the TWIX trees
}
\description{
Bootstrap samples of the Greedy-TWIX-trees.
}
\usage{
bootTWIX(formula, data=NULL,test.data=0,N=1,topN=1,subset=NULL,
                    method="deviance",topn.method="complete",
                    cluster=NULL,minsplit=30,minbucket=round(minsplit/3),
                    Devmin=0.1,level=20,score=1,tol=0.15,splitf ="entropy")
}
\arguments{
\item{formula}{formula of the form \code{y ~ x1 + x2 + \ldots},
        where \code{y} must be a factor and \code{x1,x2,\ldots} are numeric.}
\item{data}{an optional data frame containing the variables
        in the model (training data).}
\item{test.data}{a data frame containing new data.}
\item{N}{an integer giving the number of bootstrap replications.}
\item{topN}{integer vector. How many splits will be selected and at which
        level? If length 1, the same size of splits will be selected at each level.
        If length > 1, for example \code{topN=c(3,2)}, 3 splits will be chosen
        at first level, 2 splits at second level and for all next levels 1 split.}
\item{subset}{an optional vector specifying a subset of
        observations to be used.}
\item{method}{Which split points will be used? This can be \code{"deviance"}
    (default), \code{"grid"} or \code{"local"}. If the \code{method} is set to:\cr
    \code{"local"} - the program uses the local maxima of the split function(entropy),\cr
    \code{"deviance"} - all values of the entropy,\cr
    \code{"grid"} - grid points.}
\item{topn.method}{one of \code{"complete"}(default) or \code{"single"}.
    A specification of the consideration of the split points.
    If set to \code{"complete"} it uses split points from all variables,
    else it uses split points per variable.}
\item{cluster}{name of the cluster, if parallel computing will be used.}
\item{minsplit}{the minimum number of observations that must exist in a node.}
\item{minbucket}{the minimum number of observations in any terminal <leaf>
node.}
\item{Devmin}{the minimum improvement on entropy by splitting.}
\item{level}{maximum depth of the trees. If \code{level} set to 1, trees consist
        of root node.}
\item{score}{a parameter, which can be \code{1}(default) or \code{2}.
     If it is \code{2} the \emph{sort}-function will be used,\cr
     if it set to \code{1} \emph{weigth}-function will be used\cr
     \code{score = 0.25*scale(dev.tr)+0.6*scale(fit.tr)+0.15*(tree.structure)}\cr
  }
\item{tol}{parameter, which will be used, if \code{topn.method} is set to
    \code{"single"}.}
\item{splitf}{kind of the splitting function to be used. It can be one of \code{"entropy"}(default) or \code{"p-adj"}.
	If set to \code{"p-adj"}, the p-value adjusted classification tree will be build.}
}
\value{a list with the following components :
\item{call}{the call generating the object.}
\item{trees}{a list of all constructed trees, which include ID, Dev ... for each
tree.}
}
\seealso{
\code{\link{get.tree}},
\code{\link{predict.TWIX}},
\code{\link{deviance.TWIX}},\code{\link{bagg.TWIX}},
}
\examples{
data(olives)
#Tree <- bootTWIX(Region~.,data=olives,N=5)
#Tree$trees
}
\keyword{tree}
