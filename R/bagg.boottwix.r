bagg <- function(object, ...) UseMethod("bagg")

bagg.bootTWIX <- bagg.TWIX <-
    function(object, ...) bagg.default(object, ...)

bagg.default <- function(object,data=NULL,sq=1:10,...) {
    if(is.null(data))
        stop("Data not supplied")
    if (any(is.na(data)))
        stop("missing values in data")
    data <- model.frame(delete.response(terms(object$formula)),
        na.action=na.omit,data)
    if(!inherits(object, "TWIX") & !inherits(object, "bootTWIX"))
        stop("Not legitimate object")
    tpred <- predict(object,data,sq)
    apply(tpred,1,function(x) { names(sort(table(x),decreasing=TRUE)[1])})
}
