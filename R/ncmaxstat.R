##
##
## This is modified Code from maxstat Package
##
##

ncmaxstat <- function(y, x=NULL, pmethod="Lau94", minprop = 0.1, 
          maxprop=0.9, alpha = NULL, MM=0)
{
	YY<-y
	XX<-x
	N <- length(y)
	y <- y[order(x)]
	x <- my.sort(x)
	m <- MM
	
	ss <- sum(y)
	E <- m/N*ss
	
	V <- m*(N-m)/(N^2*(N-1))*(N*sum(y^2) - ss^2)
	
	Test <- abs((cumsum(y)[m] - E)/sqrt(V))
		
	STATISTIC <- max(Test)
	

	if(length(Test == STATISTIC) > 0 & length(which(Test == STATISTIC)) > 0){
		ESTIMATOR <- x[m[min(which(Test == STATISTIC))]]
	}
	else
		ESTIMATOR <- NA

#	RR <- .Call("maxstat",as.numeric(YY),as.numeric(XX),as.numeric(minprop),as.numeric(maxprop),0.0,MM,PACKAGE="TWIX")

	if(pmethod == "Lau94") {
		PVAL <- .Call("pLausen94",
				as.numeric(STATISTIC),
				as.numeric(N),
				as.numeric(minprop),
				as.numeric(maxprop),
				m,
				PACKAGE="TWIX")
#		print(PVAL)
#		cat(" STATISTIC:",STATISTIC," N: ",N," minprop: ",minprop," maxprop: ",maxprop,"\n")
#		print(m)
#		PVAL <- pLausen94(STATISTIC, N, minprop, maxprop, m=m)
#		if(!is.null(alpha))
#			QUANT <- qLausen94(alpha, N, minprop, maxprop, m=m)
	}
	RVAL <- list(statistic = STATISTIC, p.value = PVAL,
				estimate = ESTIMATOR,
				stats = Test, cuts = x[m])
	RVAL
}


pLausen94 <- function(b, N, minprop=0.1, maxprop=0.9, m=NULL){
	if(is.null(m)){
		m <- floor(N*minprop):floor(N*maxprop)
	}
	if(length(m) < 2){
			m1 <- m2 <- m
	}
	else{
		m1 <- m[1:(length(m)-1)]
		m2 <- m[2:length(m)]
	}
	t <- sqrt(1 - m1*(N-m2)/((N-m1)*m2))
	D <- sum(1/pi*exp(-b^2/2)*(t - (b^2/4 -1)*(t^3)/6))
	1 - (pnorm(b) - pnorm(-b)) + D
}

qLausen94 <- function(p, N, minprop=0.1, maxprop=0.9, m=NULL){
	test <- function(x)
	abs(pLausen94(x, N, minprop, maxprop, m) - p)
	return(optimize(test, interval=c(0,10))$minimum)
}
