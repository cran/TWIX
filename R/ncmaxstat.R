##
##
## This is slightly modified Code from maxstat Package
##
##

ncmaxstat <- function(y, x=NULL, pmethod=c("none", "Lau92", "Lau94",
          "exactGauss", "HL", "condMC", "min"), minprop = 0.1, 
          maxprop=0.9, alpha = NULL, ...)
{
  pmethod <- match.arg(pmethod)

  if (is.null(y) || is.null(x)) stop("no data given")

  if (!is.numeric(x)) {
    if (is.factor(x)) {
      if (!(is.ordered(x) || nlevels(x) == 2)) {
        warning("cannot order in x, returning NA")
        return(NA)
      }
    } else {
      warning("cannot order in x, returning NA")
      return(NA)
    }
  }

  xname <- deparse(substitute(x))
  yname <- deparse(substitute(y))  
  if (length(xname) == 1 & length(yname) == 1)
    DNAME <- paste(xname, "and", yname)
  else
    DNAME <- "y by x" 

  N <- length(y)

  y <- y[order(x)]
  x <- sort(x)
  ties <- duplicated(x)

  m <- which(!ties) - 1 
  if (minprop == 0 & maxprop==1) m <- m[2:(length(m)-1)] else {
    if (all(m < floor(N*minprop))) stop("minprop too large")
    if (all(m > floor(N*maxprop))) stop("maxprop too small")
    m <- m[m >= floor(N*minprop)]
    m <- m[m <= floor(N*maxprop)]
  }

  if(length(m) < 1) stop("no data between minprop, maxprop")

  ss <- sum(y)
  E <- m/N*ss
  V <- m*(N-m)/(N^2*(N-1))*(N*sum(y^2) - ss^2)

  Test <- abs((cumsum(y)[m] - E)/sqrt(V))

  STATISTIC <- max(Test)
  if(length(Test == STATISTIC) > 0 & length(which(Test == STATISTIC)) > 0){
  ESTIMATOR <- x[m[min(which(Test == STATISTIC))]]
  #names(STATISTIC) <- "M"
  #names(ESTIMATOR) <- c("estimated cutpoint")
  }
  else
    ESTIMATOR <- NA

  if (is.null(alpha)) QUANT <- NA

  if (pmethod == "none") {
    PVAL <- NA
    QUANT <- NA
  }
  #if (pmethod == "Lau92") {
  #  PVAL <- pLausen92(STATISTIC, minprop, maxprop)
  #  if (!is.null(alpha))
  #    QUANT <- qLausen92(alpha, minprop, maxprop)
  #}
  if (pmethod == "Lau94") {
    PVAL <- pLausen94(STATISTIC, N, minprop, maxprop, m=m)
    if (!is.null(alpha))
       QUANT <- qLausen94(alpha, N, minprop, maxprop, m=m)
  }
  #if (pmethod == "exactGauss") {
  #  PVAL <- pexactgauss(STATISTIC, x, minprop, maxprop, ...)
  #  if (!is.null(alpha))
  #     QUANT <- qexactgauss(alpha, x, minprop, maxprop, ...)
  #}
  #if (pmethod == "HL") {
  #  PVAL <- pmaxstat(STATISTIC, y, m)
  #  if (!is.null(alpha))
  #     QUANT <- qmaxstat(alpha, y, m)
  #}
  #if (pmethod == "condMC") {
  #  if (!is.null(alpha)) {
  #     maxdens <- qmaxperm(alpha, y, m, E, V, ...)
  #     QUANT <- maxdens$quant
  #     PVAL <- sum(maxdens$exdens$Prob[maxdens$exdens$T > STATISTIC])
  #  } else { 
  #    PVAL <- pmaxperm(STATISTIC, y, m, E, V, ...)
  #  }
  #}

  #if (pmethod == "min") {
  #  PVAL <- min(pLausen92(STATISTIC, minprop, maxprop), 
  #              pLausen94(STATISTIC, N, minprop, maxprop, m=m),  
  #              pexactgauss(STATISTIC, x, minprop, maxprop, ...),
  #              pmaxstat(STATISTIC, y, m))
  #  if (!is.null(alpha))
  #     QUANT <- NA
  #}
  RVAL <- list(statistic = STATISTIC, p.value = PVAL,
               method = pmethod,
               estimate = ESTIMATOR, data.name = DNAME,
               stats = Test, cuts = x[m], quant = QUANT)
  class(RVAL) <- "maxtest"
  RVAL
}


pLausen94 <- function(b, N, minprop=0.1, maxprop=0.9, m=NULL)
{
  if(is.null(m))
    m <- floor(N*minprop):floor(N*maxprop)
  m1 <- m[1:(length(m)-1)]
  m2 <- m[2:length(m)]
  t <- sqrt(1 - m1*(N-m2)/((N-m1)*m2))
  D <- sum(1/pi*exp(-b^2/2)*(t - (b^2/4 -1)*(t^3)/6))
  1 - (pnorm(b) - pnorm(-b)) + D
}

qLausen94 <- function(p, N, minprop=0.1, maxprop=0.9, m=NULL)
{
  test <- function(x)
    abs(pLausen94(x, N, minprop, maxprop, m) - p)

  return(optimize(test, interval=c(0,10))$minimum)
}
