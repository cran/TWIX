\name{sp.slave}
\alias{sp.slave}
\title{
Internal functions of TWIX.
}
\description{Trees with extra splits}
\usage{
sp.slave(rsp, m, test.data, Dmin = 0.01, minsplit = 20,
    minbucket=round(minsplit/3),topN = 1, method = "deviance",
    topn.method = "complete", level = 3, lev = 0, st = 1,
    tol = 0.1, K = 0 , oldspvar=0)
}
\arguments{
  \item{rsp}{a response variable.}
  \item{m}{a dataframe (training data).}
  \item{test.data}{a data frame containing new data.}
  \item{method}{Which split points will be used? This can be \code{"deviance"}
    (default), \code{"grid"} or \code{"local"}. If the \code{method} is set to:
    \emph{local} the program uses the local maxima of the split function(entropy),
    \emph{deviance} all values of the entropy, \emph{grid} grid points.}
  \item{topn.method}{one of \code{"complete"}(default) or \code{"single"}.
    A specification of the consideration of the split points.
    If set to \code{"complete"} it uses split points from all variables,
    else it uses split points per variable.}
  \item{minsplit}{the minimum number of observations that must exist in a
    node.}
  \item{minbucket}{the minimum number of observations in any terminal
    <leaf> node.}
  \item{Dmin}{the minimum improvement on entropy by splitting.}
  \item{topN}{a integer vector. How many splits will be selected and at which
    level? If length 1, the same size of splits will be selected at each level.
    If length > 1, for example \code{topN=c(3,2)}, 3 splits will be chosen
    at first level, 2 splits at second level and for all next levels 1 split.}
  \item{level}{maximum depth of the trees. If \code{level} set to 1, trees
    consist of root node.}
  \item{st}{step parameter for method \code{"grid"}.}
  \item{tol}{parameter, which will be used, if \code{topn.method} is set to
    \code{"single"}.}
  \item{lev}{Internal parameter.}
  \item{K}{k-fold cross-validation.}
  \item{oldspvar}{internal parameter.}
}
\seealso{
  \code{\link{TWIX}}, \code{\link{predict.TWIX}},
}
\keyword{tree}
