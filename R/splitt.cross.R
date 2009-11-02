splitt_cross <- function(sv, rsp, meth="deviance", topn=1, topn.meth="complete",
                lstep=1, test=FALSE, level=0, K=0, minbuck=1)
{
    if(is.factor(sv)) {
		split_cat <- .Call("split_cat",
					rsp, sv,
					as.integer(topn+1),test,
					as.numeric(K), as.numeric(minbuck),
					PACKAGE="TWIX")
		split_cat
    }
	else if( meth != 2){
	#####
	##  cross-validation
			n <- length(rsp)
            xval <- trunc(n/(n*K))
            xgr <- 1:xval
            s <- sample(rep(xgr,length=n),n)
            NN <- n-table(s)
            split_end <- .Call("split_cross",
                        as.numeric(sv),
                        as.integer(rsp),
                        as.integer(NN),
                        as.integer(fullrks(sv)),
                        as.integer(s),
                        as.integer(xval),
						as.integer(minbuck),
                        PACKAGE="TWIX")
        if (meth == 1) {
            if (length(split_end[[1]]) > 3) {
                dd <- ww <- vector()
                m <- length(split_end[[1]])
                id.l <- 1+
                    which((split_end[[1]][2:(m-1)] > split_end[[1]][1:(m-2)]) &
                        (split_end[[1]][2:(m-1)] > split_end[[1]][3:m]))
                if(sum(split_end[[1]][1] > split_end[[1]]) >= m)
                    id.l <- c(id.l,2)
                if(sum(split_end[[1]][m] > split_end[[1]]) >= m)
                    id.l <- c(id.l,(m+1))
                dd <- split_end[[1]][id.l]
                sc <- split_end[[3]][id.l]
                ww <- split_end[[2]][id.l]
                dev <- dd;wh<-ww
            }
            else {
                dev <- split_end[[1]]
				wh <- split_end[[2]]
				sc <- split_end[[3]]
                }
        }
        else if (meth == 0 && lstep > 1) {
            m <- length(split_end[[1]])
            if (m != 0 ) {
                for (i in 1:3) {
                    split_end[[i]] <- split_end[[i]][seq(1,m,lstep)]
                }
            }
            dev <- split_end[[1]]
			wh <- split_end[[2]]
			sc <- split_end[[3]]
        }
        else {
            dev <- split_end[[1]]
            wh <- split_end[[2]]
            sc <- split_end[[3]]
            }
        if(test || length(dev) <= topn){
            list(dev=dev,globD=split_end[[4]],which=wh,score=sc)
        }
        else {
            score<-0.7*dev/max(dev) + 0.3*sc/max(sc)
            id <- my.sort(score,decreasing=TRUE,index.return = TRUE)$ix
            if(topn != 0){
                dev <- na.omit(dev[id[1:(topn+1)]])
                wh <- na.omit(wh[id[1:(topn+1)]])
                sc <- na.omit(sc[id[1:(topn+1)]])
            }
            else{
                dev <- na.omit(dev[id])
                wh <- na.omit(wh[id])
                sc <- na.omit(sc[id])
            }
            list(dev=dev,globD=split_end[[4]],which=wh,score=sc)
            }
        }
	else if(meth == 2){
#####
##  grid-method
		n <- length(rsp)
		d <- rep(1,n)
		wh <- seq(min(sv),max(sv),length.out=n)
		if(topn > n)
			topn <- n
		id <- seq(0,n,by=n/(topn+1))
		n_id <- length(id)
		if(topn < n) {
			dev <- d[id[2:(n_id-1)]]
			wh <- wh[id[2:(n_id-1)]]
		}
		else {
			dev <- d[id[2:n_id]]
			wh <- wh[id[2:n_id]]
		}
		list(dev=dev,globD=10,which=wh)
	}
}
