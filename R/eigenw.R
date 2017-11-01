# Copyright 2002-12 by Roger Bivand 
#

eigenw <- function(listw, quiet=NULL)
{
	if(!inherits(listw, "listw")) stop("not a listw object")
        if (is.null(quiet)) quiet <- !get("verbose", envir = .spdepOptions)
        stopifnot(is.logical(quiet))

	w <- listw2mat(listw)
        sym <- all(w == t(w))
	e <- eigen(w, symmetric=sym, only.values=TRUE)$values
	if (!quiet) {
		cat("Largest eigenvalue:", 
# modified 110414 RSB
		if(is.complex(e)) {max(Re(e[which(Im(e) == 0)]))} else max(e),
		"Sum of eigenvalues:", sum(e), "\n")
	}
	e
}

griffith_sone <- function(P, Q, type="rook") {
    stopifnot(P >= 1)
    stopifnot(Q >= 1)
    p <- seq(1:P)
    q <- seq(1:Q)
    if (type=="rook") {
        res0 <- outer((2*cos((pi*p)/(P+1))), (2*cos((pi*q)/(Q+1))), FUN="+")
    } else {
        e2a <- outer((cos((pi*p)/(P+1))), (cos((pi*q)/(Q+1))), FUN="+")
        e2b <- outer((2*cos((pi*p)/(P+1))), (cos((pi*q)/(Q+1))), FUN="*")
        res <- 2*(e2a+e2b)
    }
    res <- sort(c(res0), decreasing=TRUE)
    res
}

subgraph_eigenw <- function(nb, glist=NULL, style="W", zero.policy=NULL,
    quiet=NULL) {
    if(!inherits(nb, "nb")) stop("Not a neighbours list")
    if (is.null(quiet)) quiet <- !get("verbose", envir = .spdepOptions)
    stopifnot(is.logical(quiet))
    if (is.null(zero.policy))
        zero.policy <- get("zeroPolicy", envir = .spdepOptions)
    stopifnot(is.logical(zero.policy))
    if (!(style %in% c("W", "B", "C", "S", "U", "minmax")))
        stop(paste("Style", style, "invalid"))
    can.sim <- FALSE
    if (style %in% c("W", "S"))
        can.sim <- can.be.simmed(nb2listw(nb, glist=glist, style=style))
    nc <- n.comp.nb(nb)
    t0 <- table(nc$comp.id)
    elist <- vector(mode="list", length=length(t0))
    singleton <- names(t0)[which(t0 == 1)]
    if (length(singleton) > 0) elist[as.integer(singleton)] <- 0.0
    doubles <- names(t0)[which(t0 == 2)]
    if (length(doubles) > 0) {
        for (i in doubles) elist[[as.integer(i)]] <- c(1.0, -1.0)
    }
    rest <- which(sapply(elist, is.null))
    for (i in rest) {
        nbi <- subset(nb, nc$comp.id == i)
        gli <- NULL
        if (!is.null(glist)) gli <- subset(glist, nc$comp.id == i)
        if (can.sim) {
            elist[[i]] <- eigenw(similar.listw(nb2listw(nbi, glist=gli,
                style=style)))
        } else {
            elist[[i]] <- eigenw(nb2listw(nbi, glist=gli, style=style))
        }
    }
    eout <- sort(unlist(elist))
    if (length(eout) != length(nb))
        stop("length mismatch, eout:", length(eout))
    eout
}
