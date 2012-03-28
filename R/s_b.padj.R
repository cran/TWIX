b <- function(n,NF,cp,Svar,Sp,d,d.t,ftr,ftest,node.cl,yprob,Obsn,dis,kst,ddD,L,R){
	t <- 1
	NR <- length(R)
	NL <- length(L)
	tree <- vector(mode="list",NL*NR)
	for (i in 1:NL) {
		for (j in 1:NR) {
			CP <- (1/NF)*(ddD - sum(L[[i]]$risk[L[[i]]$id == 0]) - sum(R[[j]]$risk[R[[j]]$id == 0]))/(sum(c(L[[i]]$id == 0,R[[j]]$id == 0))-1)
			if((length(L[[i]]$id) > 1 | length(R[[j]]$id) > 1) && CP >= cp){
				tree[[t]] <- 
					list(id = c(n,L[[i]]$id,R[[j]]$id),
						dev = d+L[[i]]$dev+R[[j]]$dev,
						dev.test = L[[i]]$dev.test+R[[j]]$dev.test,
						fit.tr = ftr+L[[i]]$fit.tr+R[[j]]$fit.tr,
						fit.test = ftest+L[[i]]$fit.test+R[[j]]$fit.test,
						Splitvar = c(Svar,L[[i]]$Splitvar,R[[j]]$Splitvar),
						Splitp = c(Sp[1],L[[i]]$Splitp,R[[j]]$Splitp),
						Pred.class=node.cl,
						y.prob = c(yprob,L[[i]]$y.prob,R[[j]]$y.prob),
						Obs = c(Obsn,L[[i]]$Obs,R[[j]]$Obs),
						dist = c(dis,L[[i]]$dist,R[[j]]$dist),
						ks.t = c(kst,L[[i]]$ks.t,R[[j]]$ks.t),
						risk = c(ddD,L[[i]]$risk,R[[j]]$risk))
			}
			else if((L[[i]]$Pred.class == R[[j]]$Pred.class) && CP <= cp) {
				tree[[t]] <- 
						list(id = 0, dev = 0,
							dev.test = d.t,
							fit.tr = ftr+L[[i]]$fit.tr+R[[j]]$fit.tr,
							fit.test = ftest+L[[i]]$fit.test+R[[j]]$fit.test,
							Splitvar = Svar, Splitp=Sp[1],
							Pred.class=node.cl, y.prob=yprob,
							Obs = sum(c(L[[i]]$Obs,R[[j]]$Obs)),
							dist = 0, ks.t = 0, risk=ddD)
			}
			else {
				tree[[t]] <- 
						list(id = c(n,L[[i]]$id,R[[j]]$id),
							dev = d+L[[i]]$dev+R[[j]]$dev,
							dev.test = L[[i]]$dev.test+R[[j]]$dev.test,
							fit.tr = ftr+L[[i]]$fit.tr+R[[j]]$fit.tr,
							fit.test = ftest+L[[i]]$fit.test+R[[j]]$fit.test,
							Splitvar = c(Svar,L[[i]]$Splitvar,R[[j]]$Splitvar),
							Splitp = c(Sp[1],L[[i]]$Splitp,R[[j]]$Splitp),
							Pred.class=node.cl,
							y.prob = c(yprob,L[[i]]$y.prob,R[[j]]$y.prob),
							Obs = c(Obsn,L[[i]]$Obs,R[[j]]$Obs),
							dist = c(dis,L[[i]]$dist,R[[j]]$dist),
							ks.t = c(kst,L[[i]]$ks.t,R[[j]]$ks.t),
							risk = c(ddD,L[[i]]$risk,R[[j]]$risk))
			}
			t<-t+1
		}
	}
	tree
}


s <- function(a,NF,cp) {
		id <- vector(mode="list")
        for (k in 1:length(a$split)) {
			aspl <- a$split[[k]]
			alf <- a$left[[k]]
			arf <- a$right[[k]]
            id <- c(id,
                b(k,NF,cp,aspl$Splitvar, aspl$Splitp, aspl$Dev,
                        aspl$Dev.test, 0, 0,
                        aspl$Pred.class, max(aspl$Prob),
						0, aspl$dist, aspl$ks.t, aspl$Obs-aspl$fit.tr,
                    if(length(alf) == 3)
                        s(alf,NF,cp)
                    else
                        list(list(id=0,Splitvar=0,Splitp=0,dev=0,
                        dev.test=alf$Dev.test,
                        fit.tr=alf$fit.tr, fit.test=alf$fit.test,
                        Pred.class=alf$Pred.class, y.prob=max(alf$Prob),
						Obs=alf$Obs, dist=0, ks.t=0, risk=alf$Obs-alf$fit.tr)),
                    if(length(arf) == 3)
                        s(arf,NF,cp)
                    else
                        list(list(id=0,Splitvar=0,Splitp=0,dev=0,
                        dev.test=arf$Dev.test,
                        fit.tr=arf$fit.tr, fit.test=arf$fit.test,
                        Pred.class=arf$Pred.class, y.prob=max(arf$Prob),
						Obs=arf$Obs, dist=0, ks.t=0, risk=arf$Obs-arf$fit.tr))
                )
            )
		}
	id
}



b_padj <- function(n,Svar,Sp,d,d.t,ftr,ftest,node.cl,yprob,Obsn,dis,kst,ddD,L,R){
	t <- 1
	NR <- length(R)
	NL <- length(L)
	tree <- vector(mode="list",NL*NR)
	for (i in 1:NL) {
		for (j in 1:NR) {
			if(length(L[[i]]$id) > 1 | length(R[[j]]$id) > 1){             
				tree[[t]] <- 
					list(id = c(n,L[[i]]$id,R[[j]]$id),
						dev = L[[i]]$dev+R[[j]]$dev,
						dev.test = L[[i]]$dev.test+R[[j]]$dev.test,
						fit.tr = ftr+L[[i]]$fit.tr+R[[j]]$fit.tr,
						fit.test = ftest+L[[i]]$fit.test+R[[j]]$fit.test,
						Splitvar = c(Svar,L[[i]]$Splitvar,R[[j]]$Splitvar),
						Splitp = c(Sp[1],L[[i]]$Splitp,R[[j]]$Splitp),
						Pred.class=node.cl,
						y.prob = c(yprob,L[[i]]$y.prob,R[[j]]$y.prob),
						Obs = c(Obsn,L[[i]]$Obs,R[[j]]$Obs),
						dist = c(dis,L[[i]]$dist,R[[j]]$dist),
						ks.t = c(kst,L[[i]]$ks.t,R[[j]]$ks.t),
						dD = c(ddD,L[[i]]$dD,R[[j]]$dD))
			}
			else if(L[[i]]$Pred.class == R[[j]]$Pred.class){
				tree[[t]] <- 
					list(id = 0, dev = d, dev.test = d.t,
						fit.tr = ftr+L[[i]]$fit.tr+R[[j]]$fit.tr,
						fit.test = ftest+L[[i]]$fit.test+R[[j]]$fit.test,
						Splitvar = Svar, Splitp = Sp[1],
						Pred.class=node.cl, y.prob=yprob,
						Obs = L[[i]]$Obs+R[[j]]$Obs,
						dist = 0, ks.t = 0, dD = 0)
			}
			else {
				tree[[t]] <- 
					list(id = c(n,L[[i]]$id,R[[j]]$id),
						dev = L[[i]]$dev+R[[j]]$dev,
						dev.test = L[[i]]$dev.test+R[[j]]$dev.test,
						fit.tr = ftr+L[[i]]$fit.tr+R[[j]]$fit.tr,
						fit.test = ftest+L[[i]]$fit.test+R[[j]]$fit.test,
						Splitvar = c(Svar,L[[i]]$Splitvar,R[[j]]$Splitvar),
						Splitp = c(Sp[1],L[[i]]$Splitp,R[[j]]$Splitp),
						Pred.class = node.cl,
						y.prob = c(yprob,L[[i]]$y.prob,R[[j]]$y.prob),
						Obs = c(Obsn,L[[i]]$Obs,R[[j]]$Obs),
						dist = c(dis,L[[i]]$dist,R[[j]]$dist),
						ks.t = c(kst,L[[i]]$ks.t,R[[j]]$ks.t),
						dD = c(ddD,L[[i]]$dD,R[[j]]$dD))
			}
			t <- t+1
		}
	}
    tree
}


s_padj <- function(a) {
	id <- list()
	for (k in 1:length(a$split)) {
		aspl <- a$split[[k]]
		alf <- a$left[[k]]
		arf <- a$right[[k]]
		id <- c(id,
			b_padj(k, aspl$Splitvar, aspl$Splitp, aspl$dist,
					aspl$Dev.test, 0, 0, aspl$Pred.class,
					max(aspl$Prob), 0, aspl$dist, aspl$ks.t, aspl$dD,
				if(length(alf) == 3)
					s_padj(alf)
				else
					list(list(id=0,Splitvar=0,Splitp=0,
						dev=alf$Dev,
						dev.test=alf$Dev.test,
						fit.tr=alf$fit.tr, fit.test=alf$fit.test,
						Pred.class=alf$Pred.class,
						y.prob=max(alf$Prob), Obs=alf$Obs,
						dist=0, ks.t=0, dD=0)),
				if(length(arf) == 3)
					s_padj(arf)
				else
					list(list(id=0,Splitvar=0,Splitp=0,
						dev=arf$Dev,
						dev.test=arf$Dev.test,
						fit.tr=arf$fit.tr, fit.test=arf$fit.test,
						Pred.class=arf$Pred.class, y.prob=max(arf$Prob),
						Obs=arf$Obs, dist=0, ks.t=0, dD=0))
			)
		)
	}
	id
}


#####
##### optimized Versions for bagging
#####

b_bag <- function(n,Svar,Sp,d,ftr,node.cl,L,R){
	t <- 1
	NR <- length(R)
	NL <- length(L)
	tree <- vector(mode="list",NL*NR)
	for (i in 1:NL){
		for (j in 1:NR){
			if(length(L[[i]]$id) > 1 | length(R[[j]]$id) > 1){
				tree[[t]] <- 
						list(id = c(n,L[[i]]$id,R[[j]]$id),
							dev = d+L[[i]]$dev+R[[j]]$dev,
							fit.tr = ftr+L[[i]]$fit.tr+R[[j]]$fit.tr,
							Splitvar = c(Svar,L[[i]]$Splitvar,R[[j]]$Splitvar),
							Splitp = c(Sp[1],L[[i]]$Splitp,R[[j]]$Splitp),
							Pred.class = node.cl)
			}
			else if(L[[i]]$Pred.class == R[[j]]$Pred.class){
				tree[[t]]<-
						list(id = 0, dev = 0, 
							fit.tr = ftr+L[[i]]$fit.tr+R[[j]]$fit.tr,
							Splitvar = Svar, Splitp = Sp[1], Pred.class = node.cl)
			}
			else{
				tree[[t]] <- 
						list(id = c(n,L[[i]]$id,R[[j]]$id),
							dev = d+L[[i]]$dev+R[[j]]$dev,
							fit.tr = ftr+L[[i]]$fit.tr+R[[j]]$fit.tr,
							Splitvar = c(Svar,L[[i]]$Splitvar,R[[j]]$Splitvar),
							Splitp = c(Sp[1],L[[i]]$Splitp,R[[j]]$Splitp),
							Pred.class = node.cl)
			}
			t<-t+1
		}
	}
	tree
}

s_bag <- function(a) {
	id <- list()
	for (k in 1:length(a$split)) {
		aspl <- a$split[[k]]
		id <- c(id,
			b_bag(k,aspl$Splitvar,aspl$Splitp,aspl$Dev,0,aspl$Pred.class,
				if(length(a$left[[k]]) == 3)
					s_bag(a$left[[k]])
				else
					list(list(id=0,Splitvar=0,Splitp=0,dev=0,
					fit.tr=a$left[[k]]$fit.tr,
					Pred.class=a$left[[k]]$Pred.class)),
				if(length(a$right[[k]]) == 3)
					s_bag(a$right[[k]])
				else
					list(list(id=0,Splitvar=0,Splitp=0,dev=0,
					fit.tr=a$right[[k]]$fit.tr,
					Pred.class=a$right[[k]]$Pred.class))
			)
		)
	}
	id
}
	
b_bag_padj <- function(n, Svar, Sp, d, ftr, node.cl, L, R){
	t <- 1
	NR <- length(R)
	NL <- length(L)
	tree <- vector(mode="list",NL*NR)
	for (i in 1:NL) {
		for (j in 1:NR) {
			if(length(L[[i]]$id) > 1 | length(R[[j]]$id) > 1){             
				tree[[t]] <- 
					list(id = c(n,L[[i]]$id,R[[j]]$id),
						dev = L[[i]]$dev+R[[j]]$dev,
						fit.tr = ftr+L[[i]]$fit.tr+R[[j]]$fit.tr,
						Splitvar = c(Svar,L[[i]]$Splitvar,R[[j]]$Splitvar),
						Splitp = c(Sp[1],L[[i]]$Splitp,R[[j]]$Splitp),
						Pred.class=node.cl)
			}
			else if(L[[i]]$Pred.class == R[[j]]$Pred.class){
				tree[[t]] <- 
					list(id = 0, dev = d,
						fit.tr = ftr+L[[i]]$fit.tr+R[[j]]$fit.tr,
						Splitvar = Svar, Splitp = Sp[1], Pred.class=node.cl)
			}
			else {
				tree[[t]] <- 
					list(id = c(n,L[[i]]$id,R[[j]]$id),
						dev = L[[i]]$dev+R[[j]]$dev,
						fit.tr = ftr+L[[i]]$fit.tr+R[[j]]$fit.tr,
						Splitvar = c(Svar,L[[i]]$Splitvar,R[[j]]$Splitvar),
						Splitp = c(Sp[1],L[[i]]$Splitp,R[[j]]$Splitp),
						Pred.class = node.cl)
			}
			t<-t+1
		}
	}
    tree
}

s_bag_padj <- function(a) {
	id <- list()
	for (k in 1:length(a$split)) {
		aspl <- a$split[[k]]
		id <- c(id,
			b_bag_padj(k, aspl$Splitvar, aspl$Splitp, aspl$dist, 0, aspl$Pred.class,
				if(length(a$left[[k]]) == 3)
					s_bag_padj(a$left[[k]])
				else
					list(list(id=0,Splitvar=0,Splitp=0,
						dev=a$left[[k]]$Dev,
						fit.tr=a$left[[k]]$fit.tr,
						Pred.class=a$left[[k]]$Pred.class)),
				if(length(a$right[[k]]) == 3)
					s_bag_padj(a$right[[k]])
				else
					list(list(id=0,Splitvar=0,Splitp=0,
						dev=a$right[[k]]$Dev,
						fit.tr=a$right[[k]]$fit.tr,
						Pred.class=a$right[[k]]$Pred.class))
                )
		)
	}
	id
}




