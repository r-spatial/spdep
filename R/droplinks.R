# Copyright 2001-24 by Roger Bivand 
#

droplinks <- function(nb, drop, sym=TRUE) {
# class to inherits Jari Oksanen 080603
  	if (!inherits(nb, "nb")) stop("not a neighbours list")
	n <- length(nb)
	cnb <- card(nb)
	if (n < 1) stop("non-positive length of nb")
	if (is.logical(drop)) {
		if(length(drop) != n) stop("Argument lengths differ")
		idrop <- which(drop == TRUE)
	} else if(is.character(drop)) {
		row.names <- as.character(attr(nb, "region.id"))
		idrop <- match(drop, row.names)
		if(any(is.na(idrop))) stop("Region to drop not found")
	} else {
		idrop <- match(drop, 1:n)
		if(any(is.na(idrop))) stop("Region to drop not found")
	}
	if((attr(nb, "sym") == FALSE) && (sym == TRUE)) {
		warning("setting sym to FALSE")
		sym <- FALSE
	}
	for (i in idrop) {
		if (sym && cnb[i] > 0) {
			for (j in nb[[i]])
				nb[[j]] <- nb[[j]][nb[[j]] != i]
		}
		nb[[i]] <- 0L
	}
	nb <- sym.attr.nb(nb)
        cans <- card(nb)
        if (get.NoNeighbourOption()) {
            if (any(cans == 0L)) warning("some observations have no neighbours")
        }
        NE <- n + sum(cans)
        if (get.SubgraphOption() && get.SubgraphCeiling() > NE) {
            ncomp <- n.comp.nb(nb)
            attr(nb, "ncomp") <- ncomp
            if (ncomp$nc > 1) warning("neighbour object has ", ncomp$nc, " sub-graphs")
        }
	nb
}

addlinks1 <- function(nb, from, to, sym=TRUE) {
  	if (!inherits(nb, "nb")) stop("not a neighbours list")
        stopifnot(length(from) == 1L)
	n <- length(nb)
	cnb <- card(nb)
	if (n < 1) stop("non-positive length of nb")
	row.names <- as.character(attr(nb, "region.id"))
	if (is.character(from)) {
		ifrom <- match(from, row.names)
		if(any(is.na(ifrom))) stop("from-region not found")
	} else {
		ifrom <- match(from, 1:n)
		if (any(is.na(ifrom))) stop("from-region not found")
	}
	if (is.character(to)) {
		ito <- match(to, row.names)
		if (any(is.na(ito))) stop("to-region not found")
	} else {
		ito <- match(to, 1:n)
		if(any(is.na(ito))) stop("to-region drop not found")
	}
	if ((attr(nb, "sym") == FALSE) && (sym == TRUE)) {
		warning("setting sym to FALSE")
		sym <- FALSE
	}
        orig <- nb[[ifrom]]
        orig <- orig[orig > 0L]
        nb[[ifrom]] <- as.integer(sort(unique(c(orig, ito))))
        if (sym) {
		for (i in ito) {
			orig <- nb[[i]]
			orig <- orig[orig > 0L]
			nb[[i]] <- as.integer(sort(unique(c(orig, ifrom))))
		}
        }
	nb <- sym.attr.nb(nb)
        NE <- n + sum(card(nb))
        if (get.SubgraphOption() && get.SubgraphCeiling() > NE) {
            ncomp <- n.comp.nb(nb)
            attr(nb, "ncomp") <- ncomp
            if (ncomp$nc > 1) warning("neighbour object has ", ncomp$nc, " sub-graphs")
        }
	nb
}
