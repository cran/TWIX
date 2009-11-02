

get.tree <- function (m.tr, n=1, id=NULL) {
    I <- 0
    if(attr(m.tr,"class")[1] == "TWIX"){
        m.tree <- m.tr[[5]]
		if(is.null(id))
			tree.id <- m.tr[[6]][[n]][[1]]
		else
			tree.id <- id
    }
	else{
		if (attr(m.tr,"class")[1] == "bootTWIX" || attr(m.tr,"class")[1] == "bundlTWIX")
			m.tree <- m.tr[[3]][[n]]
		if(is.null(id))
			tree.id <- m.tr[[4]][[n]][[1]]
		else
			tree.id <- id
	}
    if(tree.id[1] != 0){
        tree <- list(.Call("get_tree",m.tree,as.integer(tree.id),as.integer(I),PACKAGE="TWIX"),
						formula=m.tr[[1]])
    }
    else{
        Obs <- m.tree$split[[1]]$Obs + length(m.tr$Bad.id)
        prob <- m.tree$split[[1]]$Prob
        yval <- m.tree$split[[1]]$Pred.class
        tree <- list(id=0,formula=m.tr$formula,Obs=Obs,Dev=0,Prob=prob,Pred.class=yval)
    }
	if(attr(m.tr,"class")[2] == "p-adj"){
		attr(tree,"class") <- "padj.tree"
	}
	else{
		attr(tree,"class") <- "single.tree"
	}
	tree
}




