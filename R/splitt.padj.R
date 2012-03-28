splitt.padj <- function(sv, rsp, minprop=0.1, maxprop=0.9, minbuck=1, test=FALSE){
	err <- .Call("maxstat",
			as.numeric(sv), rsp,
			as.numeric(minprop),
			as.numeric(maxprop),
			as.logical(test),
			as.integer(minbuck),
			PACKAGE="TWIX")
	err
}