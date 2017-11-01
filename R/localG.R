# Copyright 2001-3 by Roger Bivand 
#

localG <- function(x, listw, zero.policy=NULL, spChk=NULL, return_internals=FALSE) {
	if (!inherits(listw, "listw"))
		stop(paste(deparse(substitute(listw)), "is not a listw object"))
	if (!is.numeric(x))
		stop(paste(deparse(substitute(x)), "is not a numeric vector"))
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
        stopifnot(is.vector(x))
	if (any(is.na(x))) stop(paste("NA in ", deparse(substitute(x))))
	n <- length(listw$neighbours)
	if (n != length(x))stop("Different numbers of observations")
	if (is.null(spChk)) spChk <- get.spChkOption()
	if (spChk && !chkIDs(x, listw))
		stop("Check of data and weights ID integrity failed")
	gstari <- FALSE
	if (!is.null(attr(listw$neighbours, "self.included")) &&
		attr(listw$neighbours, "self.included")) gstari <- TRUE
	lx <- lag.listw(listw, x, zero.policy=zero.policy)
	if (gstari) {
		xibar <- rep(mean(x), n)
		si2 <- rep(sum(scale(x, scale=FALSE)^2)/n, n)
	} else {
		xibar <- (rep(sum(x),n) - x) / (n - 1)
		si2 <- ((rep(sum(x^2),n) - x^2) / (n - 1)) - xibar^2
	}
	Wi <- sapply(listw$weights, sum)
	S1i <- sapply(listw$weights, function(x) sum(x^2))
        EG <- Wi*xibar
	res <- (lx - EG)
	if (gstari) {
                VG <- si2*((n*S1i - Wi^2)/(n-1))
	} else {
                VG <- si2*(((n-1)*S1i - Wi^2)/(n-2))
	}
        res <- res / sqrt(VG)
        if (return_internals) {
          if (gstari) {
            attr(res, "internals") <- cbind(G=lx/sum(c(x)),
              EG=EG/sum(c(x)), VG=VG/sum(c(x))^2)
          } else {
            attr(res, "internals") <- cbind(G=lx/(sum(c(x))-c(x)),
              EG=EG/(sum(c(x))-c(x)), VG=VG/(sum(c(x))-c(x))^2)
          }
	}
        attr(res, "gstari") <- gstari
	attr(res, "call") <- match.call()
	class(res) <- "localG"
	res
}
