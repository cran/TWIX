deviance.bootTWIX  <- function(object, ...) deviance.TWIX(object, ...)

deviance.TWIX <- function(object,type="training",...) {
    typ <- charmatch(type,c("training", "test", "all"))
    if (typ == 1)
        sapply(object$trees,function(y)y$dev)
    else if (typ == 2)
        sapply(object$trees,function(y)y$dev.test)
    else if(typ == 3)
        cbind(sapply(object$trees,function(y)y$dev),
                sapply(object$trees,function(y)y$dev.test))
}

