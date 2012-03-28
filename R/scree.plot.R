scree.plot <-function(formula, data = NULL, bars = TRUE, col = "grey", type = "b", pch = 16, ylim = c(0,1),...)
{
    call <- match.call()
    scm <- match.call(expand.dots = FALSE)
    scm$xval <- scm$runs <- scm$data <- scm$minbuck <- NULL
    scm <- model.frame(formula,data)
    rsp <- model.extract(scm, "response")
    S <- lapply(scm[2:ncol(scm)],splitt,rsp)
    S <- sapply(S,function(x) x$dev[1])
    S <- S/max(S)
    out <- sort(S,decreasing=TRUE)
    if(!bars){
        plot(out,type=type,ylim=ylim,axes=FALSE,col=col,
            ylab = expression(Delta*i(t)/Delta*i[max]),
            xlab = "Variables", main = "Screeplot",...)
        box()
        axis(2)
        axis(1,1:length(out),names(out))
    }
    else {
        barplot(out,col=col, main = "Screeplot",...)
    
    }
}
