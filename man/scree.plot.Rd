\name{scree.plot}
\alias{scree.plot}
\title{Scree-plot}
\description{
 A scree plot shows the sorted maximum decrease in impurity for each variable's value.
}
\usage{
scree.plot(formula, data = NULL, bars = TRUE, col = "grey",
            type = "b", pch = 16, ylim = c(0, 1), ...)
}
\arguments{
\item{formula}{formula of the form \code{y ~ x1 + x2 + \ldots},
    where \code{y} must be a factor and \code{x1,x2,\ldots} are numeric or factor.}
\item{data}{an optional data frame containing the variables in the model.}
\item{bars}{the type of plot: barplot or lines.}
\item{col}{The colors for lines or Bar's.}
\item{type}{the color for points.}
\item{pch}{see par(pch =...).}
\item{ylim}{the y limits of the plot.}
\item{...}{other parameters to be passed through to plotting functions.}
}
\seealso{
  \code{\link{iDevplot}},\code{\link{plot.TWIX}},\code{\link{TWIX}}}
}
\examples{
data(olives)
scree.plot(Region~.,data=olives,bars=FALSE,col=2)
scree.plot(Region~.,data=olives,bars=TRUE)
}
\keyword{tree}
