pred.value <- function(x,tree) {
	SPLV <- x[names(x) == tree$split[[1]]$Splitvar]
	if(length(tree$split[[1]]$Splitp) > 1){
		split.p <- tree$split[[1]]$Splitp
		rule <- sum(SPLV == names(split.p)[split.p == 1])
		if(length(rule) > 1) rule <- 1
	}
	else {
		rule <- as.numeric(SPLV) <= tree$split[[1]]$Splitp
	}
	if (rule) {
		if(length(tree$left$split) > 0)
			Klas <- pred.value(x=x,tree$left)
		else
			Klas <- tree$left$Pred.class
	}
	else {
		if(length(tree$right$split) > 0)
			Klas <- pred.value(x=x,tree$right)
		else
			Klas <- tree$right$Pred.class
	}
	Klas
}
