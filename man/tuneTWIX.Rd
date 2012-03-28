\name{tune.TWIX}
\alias{tune.TWIX}
\alias{tune.cp.TWIX}
\alias{tune.rpart}
\alias{tune.cp.rpart}
\alias{tune.RF}
\title{Parameter tuning via K-fold cross-validation.}
\description{
This function tunes hyperparameters: minbuck, cp for TWIX or \code{rpart}, mtry for \code{randomForest}.
}
\usage{
tune.TWIX(formula, data = NULL, minbuck = seq(5,20,by=5), 
        maxdepth=30, Devmin=0.0, splitf="deviance", 
        method="deviance", topn.method="single", 
        topN=1, tol=0.25, cp=0.0, xval=10, runs = 1, 
        trace.plot=FALSE, score=1, predict="best", 
        cluster=NULL, seed.cluster=NULL, multicore=FALSE)
		  
tune.cp.TWIX(formula, data = NULL, cp=seq(0.01,0.1,length=10),
        minbuck=1, maxdepth=30, Devmin=0.0, splitf="deviance", 
        method="deviance", topn.method="single", 
        topN=1, tol=0, xval=10, runs = 1, trace.plot=FALSE, 
        score=1, predict="best", cluster=NULL, 
        seed.cluster=NULL, multicore=FALSE)
	
tune.rpart(formula, data = NULL, minbuck=seq(5,20,by=5), 
        parms=list(split="information"), maxdepth=30, 
        cp=0.0, xval=10, runs=1, trace.plot=FALSE, 
        cluster=NULL, seed.cluster=NULL, multicore=FALSE)

tune.cp.rpart(formula, data = NULL, cp=seq(0,0.3,0.05),
        minbuck=0, parms=list(split="information"), 
        maxdepth=30, xval=10, runs=10, trace.plot=FALSE, 
        cluster=NULL, seed.cluster=NULL, multicore=FALSE)

tune.RF(formula, data = NULL, mtry = 2:(ncol(data)-1), 
        ntree=500, replace=TRUE, xval=10, runs = 1, 
        trace.plot = FALSE, cluster=NULL, seed.cluster=NULL, 
        multicore=FALSE)
}
\arguments{
\item{formula}{formula of the form \code{y ~ x1 + x2 + \ldots},
    where \code{y} must be a factor and \code{x1,x2,\ldots} are numeric or factor.}
\item{data}{an optional data frame containing the variables in the model.}
\item{minbuck}{the sampling space for parameter minbuck.}
\item{maxdepth}{set the maximum depth of the final tree.}
\item{xval}{number of cross-validations.}
\item{runs}{number of runs.}
\item{splitf}{kind of the splitting function to be used. It can be one of \code{"deviance"}(default) or \code{"p-adj"}. If set to \code{"p-adj"}, the p-value adjusted classification tree will be performed.}
\item{Devmin}{the minimum improvement on entropy by splitting. If \code{"splitf"} set to \code{"p-adj"}, \code{"Devmin"} will be the significance level alpha.}
\item{method}{Which split points will be used? This can be \code{"deviance"}
    (default), \code{"grid"} or \code{"local"}. If the \code{method} is set to:\cr
    \code{"local"} - the program uses the local maxima of the split function (entropy),\cr
    \code{"deviance"} - all values of the entropy,\cr
    \code{"grid"} - grid points.}
\item{topn.method}{one of \code{"complete"}(default) or \code{"single"}.
    A specification of the consideration of the split points.
    If set to \code{"complete"} it uses split points from all variables,
    else it uses split points per variable.}
\item{topN}{integer vector. How many splits will be selected and at which level?}
\item{tol}{parameter, which will be used, if \code{topn.method} is set to \code{"single"}.}
\item{cp}{complexity parameter. See \code{\link{rpart.control}}.}
\item{score}{Specifies the method for model selection. See \code{\link{TWIX}}.}
\item{predict}{a string that specifies what predicion mehod will be used. Possible values are "best" for a best single-tree prediction and "ensemble" for prediction based on all trees (see \code{\link{bagg}}).}
\item{cluster}{the name of the snow cluster, if parallel computing will be used}
\item{seed.cluster}{an integer to be supplied to set.seed, or NULL not to set reproducible seeds.}
\item{multicore}{a logical value for parallel execution with package \code{multicore}.}
\item{parms}{optional parameters for the splitting function. See \code{\link{rpart}}.}
\item{mtry}{Number of variables randomly sampled as candidates at each split. See \code{randomForest}.}
\item{ntree}{Number of trees to grow. See \code{randomForest}.}
\item{replace}{a logical indicating whether sampling of observations is done with or without replacement.}
\item{trace.plot}{Should trace plot be ploted?}
}
\seealso{
  \code{\link{plot.TWIX}},  \code{\link{TWIX}}
}
\examples{
    library(mlbench)
    data(PimaIndiansDiabetes2)
    Pima <- na.omit(PimaIndiansDiabetes2)

    tune.TWIX(diabetes~.,data=Pima,minbuck=c(5,10,15,20,25),runs=5)
}
\keyword{utilities}
