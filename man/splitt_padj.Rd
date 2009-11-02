\name{splitt.padj}
\alias{splitt.padj}
\title{Internal function of TWIX.}
\description{
This function compute the split-point and the corresponding p-value.}
\usage{
splitt.padj(sv, rsp, minprop=0.1, maxprop=0.9, test=FALSE)
}
\arguments{                                                                           
\item{sv}{a numeric vector of predicted variable.}
\item{rsp}{ response variable.}
\item{minprop}{ at least \code{minprop}*100\% of the observations in the
                  first group.}
\item{maxprop}{ not more than \code{maxprop}*100\% of the observations in
                  the first group.}
\item{test}{logical: return the whole information.}
}
\references{
        Lausen, B., Hothorn, T., Bretz, F. and Schmacher, M. (2004).\cr
        Assessment of Optimally Selected Prognostic Factors.\cr
        \emph{Biometrical Journal} \bold{46}, 364-374.\cr
}
\keyword{tree}
