\name{TWIX}
\alias{TWIX}
\title{
  Trees with extra splits
}
\description{Trees with extra splits}
\usage{
TWIX(formula, data=NULL, test.data=NULL, subset=NULL,
        method="deviance", topn.method="complete", minsplit=30,
        minbucket=round(minsplit/3), topN=1, splitf="deviance",
        Devmin=0.05, tol=0.15, cp=0.01, level=30, st=1, score=1,
        k=0, cluster=NULL, cl.level=1, multicore=FALSE,
        trace=TRUE, trace.plot=FALSE, ...)
}
\arguments{
  \item{formula}{formula of the form \code{y ~ x1 + x2 + \ldots},
    where \code{y} must be a factor and \code{x1,x2,\ldots} are numeric or factor.}
  \item{data}{an optional data frame containing the variables in the model (training data).}
  \item{test.data}{This can be a data frame containing new data.}
  \item{subset}{an optional vector specifying a subset of observations to be used.}
  \item{method}{Which split points will be used? This can be \code{"deviance"}
    (default), \code{"grid"} or \code{"local"}. If the \code{method} is set to:\cr
    \code{"local"} - the program uses the local maxima of the split function (entropy),\cr
    \code{"deviance"} - all values of the entropy,\cr
    \code{"grid"} - grid points.}
  \item{topn.method}{one of \code{"complete"}(default) or \code{"single"}.
    A specification of the consideration of the split points.
    If set to \code{"complete"} it uses split points from all variables,
    else it uses split points per variable.}
  \item{minsplit}{the minimum number of observations that must exist in a node.}
  \item{minbucket}{the minimum number of observations in any terminal <leaf> node.}
  \item{topN}{integer vector. How many splits will be selected and at which
    level? If length 1, the same size of splits will be selected at each level.
    If \code{length > 1}, for example \code{topN=c(3,2)}, 3 splits will be chosen
    at first level, 2 splits at second level and for all next levels 1 split.}
  \item{splitf}{kind of the splitting function to be used. It can be one of \code{"deviance"}(default) or \code{"p-adj"}. 
    If set to \code{"p-adj"}, the p-value adjusted classification tree will be performed.}
  \item{Devmin}{the minimum improvement on entropy by splitting. 
    If \code{"splitf"} set to \code{"p-adj"}, \code{"Devmin"} will be the significance level alpha.}
  \item{tol}{parameter, which will be used, if \code{topn.method} is set to \code{"single"}.}
  \item{cp}{complexity parameter.}
  \item{level}{the maximum depth of the trees. If \code{level} set to 1, trees consist of root node.}
  \item{st}{step parameter for method \code{"grid"}.}
  \item{cluster}{the name of the snow cluster, if parallel computing will be used.}
  \item{cl.level}{an internal parameter of parallel computing.}
  \item{multicore}{a logical value for parallel execution with package \code{multicore}.}
  \item{score}{Specifies the method for model selection. This can be \code{1}(default), \code{2}  or \code{3}.\cr
	If it is \code{1} the weighted correct classification rate will be used,\cr
     	if it is \code{2} the \emph{sort}-function will be used,\cr
     	if it set to \code{3} the \emph{weigth}-function will be used\cr
     	\code{score = 0.25*scale(dev.tr)+0.6*scale(fit.tr)+0.15*(structure)}
  }
  \item{k}{k-fold cross-validation of split-function. k specify the part of observations which will be take in hold-out sample (k can be (0,0.5)).}
  \item{trace}{A logical for printing a training log.}
  \item{trace.plot}{Should trace plot be ploted?}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\value{a list with the following components :
  \item{call}{the call generating the object.}
  \item{trees}{a list of all constructed trees, which include \code{ID}, \code{Dev}, \code{Fit}, \code{Splitvar}, \code{...} for each tree.}
  \item{greedy.tree}{greedy tree}
  \item{multitree}{database}
  \item{score}{score values}
}

\details{
This implementation can't handle missing values. Therefore, cases with missing values must be removed.
For p-value adjusted classification trees, continuous and binary independent descriptors are implemented 
as predictors and a response variable must be categorical with two categories.
}
\references{
	Martin Theus, Sergej Potapov and Simon Urbanek (2006).\cr
	TWIX (Talk given at the 3rd Ensemble Workshop in Munich 2006).\cr 
	\url{http://theusrus.de/Talks/Talks/TWIX.pdf}\cr


        Lausen, B., Hothorn, T., Bretz, F. and Schmacher, M. (2004).\cr
        Assessment of Optimally Selected Prognostic Factors.\cr
        \emph{Biometrical Journal} \bold{46}, 364-374.
}
\seealso{
\code{\link{get.tree}}, \code{\link{predict.TWIX}},  
\code{\link{print.single.tree}}, \code{\link{plot.TWIX}},  
\code{\link{bootTWIX}}
}
\examples{
    data(olives)

    ### train and test data
    set.seed(123)
    i <- sample(572,150)
    ic <- setdiff(1:572,i)
    training <- olives[ic,]
    test <- olives[i,]

    ###
    ### TWIX Ensemble: 729 classification trees

    TM <- TWIX(Region~.,data=training[,1:9],topN=c(9,9),method="local")
    TM$trees
    get.tree(TM,1)
    pred <- predict(TM,newdata=test,sq=1)

    ### for correct classification rate
    predict(TM,newdata=test,sq=1:36,ccr=TRUE)$CCR

    ###
    ### the p-value adjusted classification tree

    library(mlbench)
    data(PimaIndiansDiabetes2)
    Pima <- na.omit(PimaIndiansDiabetes2)


    ### train and test data
    set.seed(1111)
    N <- nrow(Pima)
    icv <- sample(N,N/3)
    itr <- setdiff(1:N,icv)
    train <- Pima[itr,]
    test <- Pima[icv,]

    ###
    ### the p-value adjusted classification tree with alpha = 0.05
   
    TMa <- TWIX(diabetes~.,data=train,splitf="p-adj",Devmin=0.05)
    get.tree(TMa)

    ### for correct classification rate
    predict(TMa,newdata=test,ccr=TRUE)$CCR
}
\keyword{tree} \keyword{classif}
