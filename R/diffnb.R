# Copyright 2001-6 by Roger Bivand 
#


diffnb <- function(x, y, verbose=NULL) {
	if (!inherits(x, "nb")) stop("not a neighbours list")
	if (!inherits(y, "nb")) stop("not a neighbours list")
        if (is.null(verbose)) verbose <- get("verbose", envir = .spdepOptions)
        stopifnot(is.logical(verbose))
	n <- length(x)
	if (n < 1) stop("non-positive length of x")
	if(n != length(y)) stop("lengths differ")
	if (any(attr(x, "region.id") != attr(y, "region.id")))
		warning("region.id differ; using ids of first list")
	ids <- attr(x, "region.id")
	res <- vector(mode="list", length=n)
	for (i in 1:n) {
		xi <- x[[i]]
		yi <- y[[i]]
		xt <- xi %in% yi
		yt <- yi %in% xi
		if (!(all(xt) && all(yt))) {
			res[[i]] <- as.integer(sort(unique(c(xi[which(!xt)],
				yi[which(!yt)]))))
			if(verbose && (res[[i]] != 0))
				cat("Neighbour difference for region id:",
				ids[i], "in relation to id:", ids[res[[i]]], "\n")
		} else res[[i]] <- 0L
	}
	class(res) <- "nb"
	attr(res, "region.id") <- attr(x, "region.id")
	attr(res, "call") <- match.call()
	res <- sym.attr.nb(res)
	res
}	
	
