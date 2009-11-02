splitt.padj <- function(sv, rsp, minprop=0.1, maxprop=0.9, test=FALSE){
	err <-.Call("maxstat",
			as.numeric(sv), rsp,
			minprop, maxprop, test,
			PACKAGE="TWIX")
	err
}