print.id.tree <- function(x,sq=1:5,...){
    if(is.null(names(x))) {
    if (length(x) < length(sq))
        sq <- 1:length(x)
        for (i in sq){
            ID <- x[[i]]$id
            DEV <- format(c(x[[i]]$dev,x[[i]]$dev.test),5)
            cat("Tree",i,"\n")
            cat(" id:",ID,"\n")
            cat(" dev:      ",DEV[1],"\n")
            cat(" dev.test: ",DEV[2],"\n")
            cat("\n")
        }
        cat("\n n =",length(x),"\n")
    } else {
        ID <- x$greedy.tree[[1]]$id
        DEV <- format(c(x$greedy.tree[[1]]$dev,x$greedy.tree[[1]]$dev.test),5)
        cat("\n")
        cat("Tree",x$id[1],"\n")
        cat(" id:",ID,"\n")
        cat(" dev:      ",DEV[1],"\n")
        cat(" dev.test: ",DEV[2],"\n")
        cat("\n")
    }
    invisible(x)
}

print.TWIX <- function(x,...){
    print(names(x)[c(2,3,5:7)])
    invisible(x)
}

print.bootTWIX <- function(x,...){
    print(names(x)[c(2,4)])
    invisible(x)
}

 
