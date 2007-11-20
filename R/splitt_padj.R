splitt_padj<-function(x,y){
	xx <- sort(x)
	ties <- duplicated(xx)
	mm <- which(!ties) - 1
	mm <- mm[mm >= floor(length(y)*0.1)]
	mm <- mm[mm <= floor(length(y)*(max(table(y)/length(y))))]
	if(length(mm) > 1){
		y <- as.numeric(y)
		err <- ncmaxstat(y,x,
				smethod="Wilcoxon",
				pmethod="Lau94",
				minprop=0.1,
				maxprop=max(table(y)/length(y)),
				alpha=0.05)
		if(length(err$p.value) < 1 || is.na(err$p.value))
			err$p.value <- 1
		Sdev<- err$p.value
		if(length(Sdev) < 1 || is.na(Sdev))
			Sdev <- 1
		Sdev<-as.numeric(Sdev)
		if(length(err$estimate) < 1 || is.na(err$estimate)){
			Swhich <- 1
		}
		else{
			Swhich <- err$estimate+
				(sort(x[as.numeric(err$estimate) < x])[1] - as.numeric(err$estimate))/2
			Swhich <- as.numeric(Swhich)
		}
		if(Sdev > 1)
			list(dev=0.0,globD=1,which=0)
		else
			list(dev=1-Sdev,globD=1,which=Swhich)
	}
	else{
		list(dev=0.0,globD=1,which=0)
	}
}