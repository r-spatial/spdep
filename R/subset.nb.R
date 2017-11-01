# Copyright 2001-8 by Roger Bivand and Yong Cai
#

subset.nb <- function(x, subset, ...) {
# class to inherits Jari Oksanen 080603
    if (!inherits(x, "nb")) stop("not a neighbours list")
    if (!is.logical(subset)) stop("subset not a logical vector")
    n <- length(x)
    if (n != length(subset))
	stop("neighours list and subset vector different lengths")
    old.ids <- 1:n
    new.ids <- match(old.ids, which(subset))
    reg.id <- subset.default(attr(x, "region.id"), subset)
    x <- sym.attr.nb(x)
    xattrs <- names(attributes(x))
    z <- subset.default(x, subset)
    nz <- length(z)
    for (i in 1:nz) {
	zi <- z[[i]]
	res <- NULL
# bug report 20050107 Yong Cai, now handles no-neighbour entities correctly
	if (!(length(zi) == 1L & zi[1] == 0)) {
	    for (j in seq(along=zi)) {
	        a <- new.ids[zi[j]]
	        if (!is.na(a)) res <- c(res, a)
	    }
	}
	if (is.null(res)) z[[i]] <- 0L
	else z[[i]] <- sort(unique(res))
    }
    attr(z, "region.id") <- reg.id
    for (i in 1:length(xattrs)) {
	if (xattrs[i] != "region.id")
	    attr(z, xattrs[i]) <- attr(x, xattrs[i])
    }
    z <- sym.attr.nb(z)
    z
}


subset.listw <- function(x, subset, zero.policy=NULL, ...) {
    if (!inherits(x, "listw")) stop("not a weights list")
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
    if (!is.logical(subset)) stop("subset not a logical vector")
    nb <- x$neighbours
    vlist <- x$weights
    if (attr(vlist, "mode") != "binary") 
	stop("Not yet able to subset general weights lists")
    style <- x$style
    n <- length(nb)
    if (n != length(subset))
	stop("neighbours list and subset vector different lengths")
    subnb <- subset.nb(x=nb, subset=subset)
    sublistw <- nb2listw(neighbours=subnb, glist=NULL, style=style,
	zero.policy=zero.policy)
    sublistw
}


