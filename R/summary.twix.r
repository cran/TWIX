summary.bootTWIX  <- function(object, ...) summary.TWIX(object, ...)

summary.TWIX <- function(object, ...) {
    training <- sapply(object$trees,function(y)y$dev)
    ccr.training <- sapply(object$trees,function(y)y$fit.tr)
    obj <- cbind(training,ccr.training)
    colnames(obj) <- c("Dev.training","CCR training")
    cat("\n")
    cat("Glob.Dev.training",object$multitree[[1]][[1]]$TD)
    cat("\n","\n")
    summary(obj)

}

