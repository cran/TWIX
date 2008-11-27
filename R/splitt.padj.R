splitt_padj <- function(sv, rsp, minprop=0.1, maxprop=0.9, test=FALSE){
	err <-.Call("maxstat",
			rsp, as.numeric(sv),
			minprop, maxprop,
			PACKAGE="TWIX")
	if(is.na(err$p.value)){
		list(dev=0.0,globD=1,which=0)
	}
	else{
		if(!test)
			list(dev=1-err$p.value,globD=1,which=err$which)
		else
			err
	}
}