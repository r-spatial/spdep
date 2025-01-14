# Copyright 2001-6, 2025 #175 by Roger Bivand 
#


diffnb <- function(x, y, verbose=NULL, legacy=TRUE) {
	if (!inherits(x, "nb")) stop("not a neighbours list")
	if (!inherits(y, "nb")) stop("not a neighbours list")
        if (is.null(verbose)) verbose <- get.VerboseOption()
        stopifnot(is.logical(verbose))
	n <- length(x)
	if (n < 1) stop("non-positive length of x")
	if(n != length(y)) stop("lengths differ")
	if (any(attr(x, "region.id") != attr(y, "region.id")))
		warning("region.id differ; using ids of first list")
	ids <- attr(x, "region.id")
	if (requireNamespace("spatialreg", quietly=TRUE) && !legacy) {
		Bx <- as(nb2listw(x, style="B", zero.policy=TRUE), 
			"CsparseMatrix")
		By <- as(nb2listw(y, style="B", zero.policy=TRUE), 
			"CsparseMatrix")
                B <- Bx != By
                res <- mat2listw(B, style="B", zero.policy=TRUE)$neighbours
	} else {
		res <- vector(mode="list", length=n)
		for (i in 1:n) {
			xi <- setdiff(x[[i]], 0L)
			yi <- setdiff(y[[i]], 0L)
#			xt <- xi %in% yi
#			yt <- yi %in% xi
#			if (!(all(xt) && all(yt))) {
#				res[[i]] <- as.integer(sort(unique(c(xi[which(!xt)],
#					yi[which(!yt)]))))
                        xy <- setdiff(xi, yi)
                        yx <- setdiff(yi, xi)
			res[[i]] <- as.integer(sort(unique(union(xy, yx))))
			if (length(res[[i]]) == 0L) res[[i]] <- 0L
			if(verbose && all(res[[i]] != 0))
				cat("Neighbour difference for region id:",
				ids[i], "in relation to id:", ids[res[[i]]], "\n")
		}
	}
	class(res) <- "nb"
	attr(res, "region.id") <- attr(x, "region.id")
	attr(res, "call") <- match.call()
	res <- sym.attr.nb(res)
        NE <- n + sum(card(res))
        if (get.SubgraphOption() && get.SubgraphCeiling() > NE) {
          ncomp <- n.comp.nb(res)
          attr(res, "ncomp") <- ncomp
          if (ncomp$nc > 1) warning("neighbour object has ", ncomp$nc, " sub-graphs")
        }
	res
}	
	
