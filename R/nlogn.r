nlogn <- function(x) { if (x == 0) 0 else x*log(x) }

fullrks<-function(m) {
	if (!is.factor(m)) .Internal(qsort(m,index.return=TRUE))[[2]] else 0
}

my.sort<-function(x,index.return=FALSE,decreasing=FALSE) {
	if(decreasing) 
		x <- -x
	y <- .Internal(qsort(x,index.return=index.return))
	if(decreasing) 
		y$x <- -y$x
	y
}