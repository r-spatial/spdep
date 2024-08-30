# Copyright 2001-8 by Roger Bivand and Yong Cai
#

subset.nb <- function(x, subset, ...) {
# class to inherits Jari Oksanen 080603
    if (!inherits(x, "nb")) stop("not a neighbours list")
    if (!is.logical(subset)) stop("subset not a logical vector")
    n <- length(x)
    input_nc <- attr(x, "ncomp")$nc
    if (is.null(input_nc) && get.SubgraphOption() && 
        get.SubgraphCeiling() > (length(x) + sum(card(x)))) {
            input_nc <- n.comp.nb(x)$nc
    }
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
	if (xattrs[i] != "region.id" && xattrs[i] != "ncomp")
	    attr(z, xattrs[i]) <- attr(x, xattrs[i])
    }
    z <- sym.attr.nb(z)
    NE <- length(z) + sum(card(z))
    if (get.SubgraphOption() && get.SubgraphCeiling() > NE) {
      ncomp <- n.comp.nb(z)
      attr(z, "ncomp") <- ncomp
      if (!is.null(input_nc) && (input_nc < ncomp$nc))
          warning("subsetting caused increase in subgraph count")
    }
    z
}


subset.listw <- function(x, subset, zero.policy=attr(x, "zero.policy"), ...) {
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
    if (!is.null(attr(x, "region.id"))) 
        attr(nb, "region.id") <- attr(x, "region.id")
    subnb <- subset.nb(x=nb, subset=subset)
    if (any(card(subnb) == 0L)) {
        if (!zero.policy) {
            warning("subsetting created no-neighbour observations, zero.policy set TRUE")
            zero.policy <- !zero.policy
        }
    }
    sublistw <- nb2listw(neighbours=subnb, glist=NULL, style=style,
	zero.policy=zero.policy)
    sublistw
}


