# Copyright 2001 by Roger Bivand 
#

moran.plot <- function(x, listw, zero.policy=NULL, spChk=NULL,
 labels=NULL, xlab=NULL, ylab=NULL, quiet=NULL, ...)
{
	if (!inherits(listw, "listw")) stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
        if (is.null(quiet)) quiet <- !get("verbose", envir = .spdepOptions)
        stopifnot(is.vector(x))
        stopifnot(is.logical(quiet))
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
	xname <- deparse(substitute(x))
	if (!is.numeric(x)) stop(paste(xname, "is not a numeric vector"))
	if (any(is.na(x))) stop("NA in X")
	n <- length(listw$neighbours)
	if (n != length(x)) stop("objects of different length")
	if (is.null(spChk)) spChk <- get.spChkOption()
	if (spChk && !chkIDs(x, listw))
		stop("Check of data and weights ID integrity failed")
	labs <- TRUE
	if (is.logical(labels) && !labels) labs <- FALSE
	if (is.null(labels) || length(labels) != n)
		labels <- as.character(attr(listw, "region.id"))
	wx <- lag.listw(listw, x, zero.policy=zero.policy)
	if (is.null(xlab)) xlab <- xname
	if (is.null(ylab)) ylab <- paste("spatially lagged", xname)
	plot(x, wx, xlab=xlab, ylab=ylab, ...)
	if (zero.policy) {
		n0 <- wx == 0.0
# bug found 100401 Paulo Grahl
                if (any(n0)) {
		    symbols(x[n0], wx[n0], inches=FALSE, 
		    circles=rep(diff(range(x))/50, length(which(n0))),
		        bg="grey", add=TRUE)
                }
	}
	xwx.lm <- lm(wx ~ x)
	abline(xwx.lm)
	abline(h=mean(wx), lty=2)
	abline(v=mean(x), lty=2)
	infl.xwx <- influence.measures(xwx.lm)
	is.inf <- which(apply(infl.xwx$is.inf, 1, any))
	points(x[is.inf], wx[is.inf], pch=9, cex=1.2)
	if (labs)
	    text(x[is.inf], wx[is.inf], labels=labels[is.inf], pos=2, cex=0.7)
	rownames(infl.xwx$infmat) <- labels
	if (!quiet) summary(infl.xwx)
	invisible(infl.xwx)
}

