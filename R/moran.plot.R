# Copyright 2001-24 by Roger Bivand 
#

moran.plot <- function(x, listw, y=NULL, zero.policy=attr(listw, "zero.policy"), spChk=NULL,
 labels=NULL, xlab=NULL, ylab=NULL, quiet=NULL, plot=TRUE, return_df=TRUE, ...)
{
	if (!inherits(listw, "listw")) stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
        if (is.null(quiet)) quiet <- !get("verbose", envir = .spdepOptions)
        stopifnot(is.vector(x))
        if (!is.null(y)) stopifnot(is.vector(y))
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
        if (!is.null(y)) {
            yname <- deparse(substitute(y))
	    if (!is.numeric(y)) stop(paste(yname, "is not a numeric vector"))
	    if (any(is.na(y))) stop("NA in Y")
	    if (n != length(y)) stop("objects of different length")
	    if (spChk && !chkIDs(y, listw))
		stop("Check of data and weights ID integrity failed")
        }
	labs <- TRUE
	if (is.logical(labels) && !labels) labs <- FALSE
	if (is.null(labels) || length(labels) != n)
		labels <- as.character(attr(listw, "region.id"))
        if (!is.null(y)) {
	    wx <- lag.listw(listw, y, zero.policy=zero.policy)
        } else {
	    wx <- lag.listw(listw, x, zero.policy=zero.policy)
        }
        if (anyNA(wx)) warning("no-neighbour observation(s) in moran.plot() - use zero.policy=TRUE")
	if (is.null(xlab)) xlab <- xname
	if (is.null(ylab)) {
            if (!is.null(y)) {
                ylab <- paste("spatially lagged", yname)
            } else {
                ylab <- paste("spatially lagged", xname)
            }
        }
	if (plot) plot(x, wx, xlab=xlab, ylab=ylab, ...)
	if (plot && zero.policy) {
		n0 <- wx == 0.0
# bug found 100401 Paulo Grahl
                if (any(n0)) {
		    symbols(x[n0], wx[n0], inches=FALSE, 
		    circles=rep(diff(range(x))/50, length(which(n0))),
		        bg="grey", add=TRUE)
                }
	}
	xwx.lm <- lm(wx ~ x)
	if (plot) abline(xwx.lm)
	if (plot) abline(h=mean(wx, na.rm=TRUE), lty=2)
	if (plot) abline(v=mean(x), lty=2)
	infl.xwx <- influence.measures(xwx.lm)
	is.inf <- apply(infl.xwx$is.inf, 1, any)
	if (plot) points(x[is.inf], wx[is.inf], pch=9, cex=1.2)
	if (plot && labs)
	    text(x[is.inf], wx[is.inf], labels=labels[is.inf], pos=2, cex=0.7)
        if (length(labels) > nrow(infl.xwx$infmat)) {
            rownames(infl.xwx$infmat) <- labels[!is.na(wx)]
        } else {
	    rownames(infl.xwx$infmat) <- labels
        }
	if (!quiet) summary(infl.xwx)
        if (return_df) {
            if (length(labels) > nrow(infl.xwx$infmat)) {
                res <- data.frame(x=x[!is.na(wx)], wx=wx[!is.na(wx)],
                    is_inf=is.inf, labels=labels[!is.na(wx)])
            } else {
                res <- data.frame(x=x, wx=wx, is_inf=is.inf, labels=labels)
            }
            res <- cbind(res, as.data.frame(infl.xwx$infmat))
            attr(res, "xname") <- xname
        } else {
            res <- infl.xwx
        }
	invisible(res)
}

