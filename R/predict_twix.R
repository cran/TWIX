predict.TWIX <- function(object,newdata,sq=1,ccr=FALSE,...) {
    if (!inherits(object, "TWIX") && !inherits(object,"single.tree"))
        stop(" object not of class TWIX or single.tree ")
    if (any(is.na(newdata)))
        stop(" missing values in newdata ")
    if (class(object) == "TWIX") {
        Dt <- model.frame(delete.response(terms(object$formula)),
            na.action=na.omit,newdata)
    }
    else {
        Dt <- model.frame(delete.response(terms(object[[2]])),
            na.action=na.omit,newdata)
    }
    ############## make 'data.frame' numeric 
    fact <- sapply(Dt, function(x) !is.null(levels(x)))
    char <- sapply(Dt, is.character)
    sDt <- 1:ncol(Dt)
    if(any(fact | char)) {
        for (j in sDt[char])
            Dt[[j]] <- as.factor(Dt[[j]])
        fact <- fact | char
        cat.levels <- lapply(Dt[fact], levels)
        for (j in sDt[fact])
            Dt[[j]] <- as.integer(Dt[[j]])
        attr(Dt, "cat.levels") <- cat.levels
    }
    ##############
    Erg <- vector()
    if (class(object) == "TWIX") {
        ldata<-split(Dt,1:nrow(Dt))
        cat.levels <- attr(Dt,"cat.levels")
        mtree <- lapply(sq,function(x,y) {get.tree(y,x)[[1]]},object)
        if(!is.list(mtree[[1]]) && mtree[[1]] == 0){
            mtree<-get.tree(object,1)
            Erg <- matrix(rep(mtree[[6]],nrow(Dt)))
        }
        else{
            Erg <- .Call("pred_TWIX",
                    ldata,
                    mtree,
                    names(Dt),
                    cat.levels,
                    length(ldata),
                    PACKAGE="TWIX")
        }
    }
    if (class(object) == "single.tree") {
        ldata<-split(Dt,1:nrow(Dt))
        cat.levels <- attr(Dt,"cat.levels")
        mtree <- list(object[[1]])
        if( mtree[[1]] == 0){
            mtree<-get.tree(object,1)
            Erg <- matrix(rep(mtree[[6]],nrow(Dt)))
        }
        else{
            Erg <- .Call("pred_TWIX",
                    ldata,
                    mtree,
                    names(Dt),
                    cat.levels,
                    length(ldata),
                    PACKAGE="TWIX")
        }
    }
    if ( ccr ) {
        n <- dim(Erg)[2]
        rate <- vector(,length=n)
        y <- newdata[names(newdata) == object$formula[[2]]]
        for (j in 1:n) {
            rate[j] <- sum(Erg[ ,j] == y[[1]])/length(y[[1]])
            }
        colnames(Erg) <- colnames(Erg, do.NULL = FALSE, prefix = "Baum.")
        Erg <- data.frame(Erg)
        Erg <- list(predicted=Erg,CCR=rate)
    }
    else {
        colnames(Erg) <- colnames(Erg, do.NULL = FALSE, prefix = "Baum.")
        Erg <- data.frame(Erg)
        }
    Erg
}
