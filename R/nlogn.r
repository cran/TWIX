nlogn <- function(x) { if (x == 0) 0 else x*log(x) }

fullrks<-function(m) {
  if (!is.factor(m)) sort(m, index.return=TRUE)$ix else 0
}

