\name{bootTWIX}
\alias{bootTWIX}
\title{
Bootstrap of the TWIX trees
}
\description{
Bootstrap samples of the Greedy-TWIX-trees or p-value adjusted TWIX-trees.
}
\usage{
bootTWIX(formula, data = NULL, nbagg = 25, topN = 1, subset = NULL, 
        comb = NULL, method = "deviance", topn.method = "complete", 
        replace = TRUE, ns=1, minsplit = 2, minbucket = round(minsplit/3), 
        splitf = "deviance", Devmin = 0.01, level = 30, tol = 0.01, 
        cluster = NULL, seed.cluster = NULL) 

}

\arguments{
\item{formula}{formula of the form \code{y ~ x1 + x2 + \ldots},
        where \code{y} must be a factor and \code{x1,x2,\ldots} are numeric.}
\item{data}{an optional data frame containing the variables in the model (training data).}
\item{nbagg}{an integer giving the number of bootstrap replications.}
\item{comb}{a list of additional model and it's prediction function for model combination, see below for some examples.}
\item{splitf}{kind of the splitting function to be used. 
	It can be one of \code{"deviance"}(default) or \code{"p-adj"}.
	If \code{splitf} set to \code{"p-adj"}, the p-value adjusted classification tree will be performed.}
\item{replace}{ Should sampling be with replacement? }
\item{ns}{ data set of size \code{ns <= nrow(data)} obtained by sampling without replacement.}
\item{Devmin}{the minimum improvement on entropy by splitting or by the p-value adjusted classification trees the significance level alpha.}
\item{topN}{integer vector. How many splits will be selected and at which
        level? If length 1, the same size of splits will be selected at each level.
        If length > 1, for example \code{topN=c(3,2)}, 3 splits will be chosen
        at first level, 2 splits at second level and for all next levels 1 split.}
\item{subset}{an optional vector specifying a subset of observations to be used.}
\item{method}{Which split points will be used? This can be \code{"deviance"}
	(default), \code{"grid"} or \code{"local"}. If the \code{method} is set to:\cr
	\code{"local"} - the program uses the local maxima of the split function(entropy),\cr
	\code{"deviance"} - all values of the entropy,\cr
	\code{"grid"} - grid points.}
\item{topn.method}{one of \code{"complete"}(default) or \code{"single"}.
	A specification of the consideration of the split points.
	If set to \code{"complete"} it uses split points from all variables,
	else it uses split points per variable.}
\item{minsplit}{the minimum number of observations that must exist in a node.}
\item{minbucket}{the minimum number of observations in any terminal <leaf> node.}
\item{level}{maximum depth of the trees. If \code{level} set to 1, trees consist of root node.}
\item{tol}{parameter, which will be used, if \code{topn.method} is set to \code{"single"}.}
\item{cluster}{the name of the cluster, if parallel computing will be used.}
\item{seed.cluster}{an integer to be supplied to set.seed, or NULL not to set reproducible seeds.}
}

\value{a list with the following components :
\item{call}{the call generating the object.}
\item{trees}{a list of all constructed trees, which include \code{ID}, \code{Dev}, \code{Fit}, \code{Splitvar}, \code{...} for each tree.}
}

\seealso{
\code{\link{TWIX}},  \code{\link{get.tree}},  
\code{\link{predict.bootTWIX}},  \code{\link{deviance.TWIX}},  
\code{\link{bagg.TWIX}}
}

\examples{
    library(ElemStatLearn)
    data(SAheart)

    ### response variable must be a factor
    SAheart$chd <- factor(SAheart$chd) 


    ### test and train data
    ###
    set.seed(1234)
    icv <- sample(nrow(SAheart),nrow(SAheart)*0.3)
    itr <- setdiff(1:nrow(SAheart),icv)
    train <- SAheart[itr,]
    test <- SAheart[icv,]


    ### Bagging with greedy decision trees as base classifier
    M1 <- bootTWIX(chd~.,data=train,nbagg=50)



    ### Bagging with the p-value adjusted classification trees (alpha is 0.01)
    ### as base classifier
    M2 <- bootTWIX(chd~.,data=train,nbagg=50,splitf="p-adj",Devmin=0.01)


    library(MASS)

    ### Double-Bagging: combine LDA and classification trees
    comb.lda <- list(model=lda, predict=function(obj, newdata)
                     predict(obj, newdata)$x)

    M3 <- bootTWIX(chd~.,data=train,nbagg=50,comb=comb.lda)


    ### Double-Bagging: combine LDA and p-value adjusted classification trees
    comb.lda <- list(model=lda, predict=function(obj, newdata)
                     predict(obj, newdata)$x)

    M4 <- bootTWIX(chd~.,data=train,nbagg=50,comb=comb.lda,
                     splitf="p-adj",Devmin=0.01)


    ### Double-Bagging: combine GLM and classification trees
    comb.glm <- list(model=function(x,...){glm(x,family=binomial,...)}, 
                     predict=function(obj, newdata) predict(obj, newdata))

    M5 <- bootTWIX(chd~.,data=train,nbagg=50,comb=comb.glm)


    ### Double-Bagging: combine GLM and p-value adjusted classification trees
    comb.glm <- list(model=function(x,...){glm(x,family=binomial,...)}, 
                     predict=function(obj, newdata) predict(obj, newdata))

    M6 <- bootTWIX(chd~.,data=train,nbagg=50,comb=comb.glm,
                     splitf="p-adj",Devmin=0.01)


    pred1 <- predict(M1,test)
    pred2 <- predict(M2,test)
    pred3 <- predict(M3,test)
    pred4 <- predict(M4,test)
    pred5 <- predict(M5,test)
    pred6 <- predict(M6,test)

    ###
    ### CCR's
  
    sum(pred1 == test$chd)/nrow(test)
    sum(pred2 == test$chd)/nrow(test)
    sum(pred3 == test$chd)/nrow(test)
    sum(pred4 == test$chd)/nrow(test)
    sum(pred5 == test$chd)/nrow(test)
    sum(pred6 == test$chd)/nrow(test)
}
\keyword{tree} \keyword{classif}
