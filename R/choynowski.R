# Copyright 2004-2010 by Roger Bivand
#

choynowski <- function(n, x, row.names=NULL, tol=.Machine$double.eps^0.5, legacy=FALSE) {
  len <- length(n)
  if (len < 1) stop("non-positive number of observations")
  res <- numeric(len)
  nsum <- sum(n)
  xsum <- sum(x)
  b <- nsum/xsum
  if (b > 1) stop("sum of cases larger than sum of populations at risk")
  E <- x*b
  type <- (n < E)
  if (legacy) {
   for (i in 1:len) {
    if(type[i]) {
      for (j in 0:n[i]) {
        xx <- (E[i]^j*exp(-E[i])) / gamma(j + 1)
        res[i] <- res[i] + xx
      }
    } else {
      xx <- 1
      x <- n[i]
      while (xx > tol) {
        xx <- (E[i]^x*exp(-E[i])) / gamma(x + 1)
        res[i] <- res[i] + xx
        x <- x + 1
      }
    }
   }
  } else {
    res <- ifelse(type, ppois(n, E), 1 - ppois(n-1, E))
  }
  if (is.null(row.names)) 
    res <- data.frame(pmap=res, type=type)
  else
    res <- data.frame(pmap=res, type=type, row.names=row.names)
  res
}
