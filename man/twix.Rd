\name{TWIX}
\alias{TWIX}
\title{
  Trees with extra splits
}
\description{Trees with extra splits}
\usage{
TWIX(formula, data = NULL, test.data = 0, subset = NULL,
        method = "deviance", topn.method = "complete", cluster = NULL,
        minsplit = 30, minbucket = round(minsplit/3), Devmin = 0.05,
	splitf = "entropy", topN = 1, level = 30, st = 1, cl.level = 2,
	tol = 0.15, score = 1,k = 0, trace.plot=FALSE,robust = FALSE, ...)
}
\arguments{
  \item{formula}{formula of the form \code{y ~ x1 + x2 + \ldots},
    where \code{y} must be a factor and \code{x1,x2,\ldots} are numeric or factor.}
  \item{data}{an optional data frame containing the variables
    in the model (training data).}
  \item{test.data}{This can be a data frame containing new data, \code{0} (default),
    or \code{"NULL"}. If set to \code{"NULL"} the bad obserations will be specified.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used.}
  \item{method}{Which split points will be used? This can be \code{"deviance"}
    (default), \code{"grid"} or \code{"local"}. If the \code{method} is set to:\cr
    \code{"local"} - the program uses the local maxima of the split function (entropy),\cr
    \code{"deviance"} - all values of the entropy,\cr
    \code{"grid"} - grid points.}
  \item{topn.method}{one of \code{"complete"}(default) or \code{"single"}.
    A specification of the consideration of the split points.
    If set to \code{"complete"} it uses split points from all variables,
    else it uses split points per variable.}
  \item{cluster}{name of the cluster, if parallel computing will be used.}
  \item{minsplit}{the minimum number of observations that must exist in a
    node.}
  \item{minbucket}{the minimum number of observations in any terminal
    <leaf> node.}
  \item{Devmin}{the minimum improvement on entropy by splitting.}
  \item{topN}{integer vector. How many splits will be selected and at which
    level? If length 1, the same size of splits will be selected at each level.
    If length > 1, for example \code{topN=c(3,2)}, 3 splits will be chosen
    at first level, 2 splits at second level and for all next levels 1 split.}
  \item{level}{maximum depth of the trees. If \code{level} set to 1, trees
    consist of root node.}
  \item{st}{step parameter for method \code{"grid"}.}
  \item{cl.level}{parameter for parallel computing.}
  \item{tol}{parameter, which will be used, if \code{topn.method} is set to
    \code{"single"}.}
  \item{score}{a parameter, which can be \code{1}(default) or \code{2}.
     If it is \code{2} the \emph{sort}-function will be used,\cr
     if it set to \code{1} \emph{weigth}-function will be used\cr
     \code{score = 0.25*scale(dev.tr)+0.6*scale(fit.tr)+0.15*(tree.structure)}\cr
  }
  \item{k}{k-fold cross-validation of split-function.
    k specify the part of observations which will be take in hold-out sample (k can be (0,0.5)).}
  \item{trace.plot}{Should trace plot be ploted?}
  \item{splitf}{kind of the splitting function to be used. It can be one of \code{"entropy"}(default) or \code{"p-adj"}.
	If set to \code{"p-adj"}, the p-value adjusted classification tree will be build.}
  \item{robust}{If set to \code{TRUE}, robust entropy estimation will be used.}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\value{a list with the following components :
  \item{call}{the call generating the object.}
  \item{trees}{a list of all constructed trees, which include ID, Dev ... 
    for each tree.}
  \item{greedy.tree}{greedy tree}
  \item{multitree}{database}
  \item{agg.id }{ vector specifying trees for aggregation.}
  \item{Bad.id}{ID-vector of bad observations from train data.}
}
\references{
		Martin Theus, Sergej Potapov and Simon Urbanek (2006).
		TWIX (Talk given at the 3rd Ensemble Workshop in Munich 2006)
		\url{http://theusrus.de/Talks/Talks/TWIX.pdf}
		
        Lausen, B., Hothorn, T., Bretz, F. and Schmacher, M. (2004).
        Assessment of Optimally Selected Prognostic Factors.
        Biometrical Journal 46, 364-374.
}
\seealso{
  \code{\link{get.tree}}, \code{\link{predict.TWIX}},
  \code{\link{print.single.tree}}, \code{\link{plot.TWIX}},
  \code{\link{deviance.TWIX}}
}
\examples{


data(olives)
set.seed(123)
i <- sample(572,150)
ic <- setdiff(1:572,i)
training <- olives[ic,]
test <- olives[i,]

TM1<-TWIX(Region~.,data=training[,1:9],topN=c(4,3),method="local")
TM1$trees
pred <- predict(TM1,newdata=test,sq=1)
predict(TM1,newdata=test,sq=1:2,ccr=TRUE)$CCR


##
## the p-value adjusted classification tree

data(iris)
TM2 <- TWIX(Species~.,data=iris,splitf="p-adj")
predict(TM2,newdata=iris,sq=1,ccr=TRUE)$CCR

}
\keyword{tree}
