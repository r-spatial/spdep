# Copyright 2001-15 by Roger Bivand 
#

spweights.constants <- function(listw, zero.policy=NULL, adjust.n=TRUE) {
    if (get.listw_is_CsparseMatrix_Option()) {
        stopifnot(is(listw, "CsparseMatrix"))
        cards <- rowSums(listw > 0)
	if (adjust.n) {
            n <- as.double(length(which(cards > 0)))
	} else {
            n <- as.double(length(cards))
        }
        c1 <- rowSums(listw)
        S0 <- sum(c1)
        S1 <- sum((listw*listw)+(listw*t(listw)))
        S2 <- sum((rowSums(listw)+colSums(listw))^2)
    } else {
	if(!inherits(listw, "listw")) stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	cards <- card(listw$neighbours)
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
	if (!zero.policy && any(cards == 0))
		stop("regions with no neighbours found")
	if (adjust.n) n <- as.double(length(which(cards > 0)))
	else n <- as.double(length(cards))
	S0 <- Szero(listw)
	S1 <- 0
	rS <- numeric(length(listw$neighbours))
	cS <- numeric(length(listw$neighbours))
	for (i in 1:length(listw$neighbours)) {
		cond <- TRUE
		if (zero.policy && cards[i] == 0) cond <- FALSE
		if (cond) {
# Luc Anselin 2006-11-11 problem with asymmetric listw
			if (cards[i] == 0)
				stop(paste("region", i,
					"has no neighbours"))
			ij <- listw$neighbours[[i]]
			wij <- listw$weights[[i]]
			rS[i] <- sum(wij)
			for (j in 1:length(ij)) {
				dij <- wij[j]
				ij.j <- ij[j]
				cS[ij.j] <- cS[ij.j] + dij
				ij.lkup <- which(listw$neighbours[[ij.j]] == i)
				if (length(ij.lkup) == 1L)
					dji <- listw$weights[[ij.j]][ij.lkup]
				else dji <- 0
				S1 <- S1 + (dij*dij) + (dij*dji)
			}
		}
	}
	S2 <- sum((rS + cS)^2)
    }
    n1 <- n - 1
    n2 <- n - 2
    n3 <- n - 3
    nn <- n*n
    list(n=n, n1=n1, n2=n2, n3=n3, nn=nn, S0=S0, S1=S1, S2=S2)
}

Szero <- function(listw) {
	sum(unlist(listw$weights))
}

lag.listw <- function(x, var, zero.policy=NULL, NAOK=FALSE, ...) {
    if (!is.logical(NAOK)) stop("NAOK must be logical")
    if (get.listw_is_CsparseMatrix_Option()) {
	if (!NAOK && any(is.na(var))) stop("NA in variable")
        stopifnot(is(x, "CsparseMatrix"))
        res <- drop(as.matrix(unname((x %*% var))))
    } else {
	listw <- x
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
	if (!inherits(listw, "listw")) stop(paste(deparse(substitute(x)),
		"is not a listw object"))
	x <- var
	if (!is.vector(c(x)) && !is.matrix(x)) stop(paste(deparse(substitute(var)),
		"not a vector or matrix"))
	if (!is.numeric(x)) stop(paste(deparse(substitute(var)),
		"not numeric"))
        storage.mode(x) <- "double"
	n <- length(listw$neighbours)
	cardnb <- card(listw$neighbours)
	if (is.null(dim(x))) {
		if (length(x) != n) stop("object lengths differ")
		res <- .Call("lagw", listw$neighbours, listw$weights,
			x, as.integer(cardnb),
			as.logical(zero.policy), naok=NAOK, PACKAGE="spdep")
	} else {
		if (nrow(x) != n) stop("object lengths differ")
		res <- matrix(0, nrow=nrow(x), ncol=ncol(x))
		for (i in 1:ncol(x)) {
			res[,i] <- .Call("lagw", listw$neighbours,
				listw$weights, x[,i],
				as.integer(cardnb), as.logical(zero.policy),
				naok=NAOK, PACKAGE="spdep")

		}
	} 
    }
    if (any(is.na(res))) warning("NAs in lagged values")
    res
}

listw2U <- function(listw) {
    if (get.listw_is_CsparseMatrix_Option()) {
        stopifnot(is(listw, "CsparseMatrix"))
        res <- (listw+t(listw))/2
    } else {
	if (!inherits(listw, "listw")) stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	nb <- listw$neighbours
	wts <- listw$weights
	style <- paste(listw$style, "U", sep="")
	sym <- is.symmetric.nb(nb, FALSE, TRUE)
	n <- length(listw$neighbours)
	cardnb <- card(listw$neighbours)
	nlist <- vector(mode="list", length=n)
	attr(nlist, "region.id") <- attr(nb, "region.id")
	class(nlist) <- "nb"
	vlist <- vector(mode="list", length=n)
	attr(vlist, as.character(style)) <- TRUE
	if (sym) {
		nlist <- vector(mode="list", length=n)
		attr(nlist, "region.id") <- attr(nb, "region.id")
		class(nlist) <- "nb"
		for (i in 1:n) {
			inb <- nb[[i]]
			nlist[[i]] <- inb
			iwt <- wts[[i]]
			icd <- cardnb[i]
			if (icd > 0) {
			    for (j in 1:icd) {
				vlist[[i]][j] <- 0.5 *
				(iwt[j]+wts[[inb[j]]][which(nb[[inb[j]]] == i)])
			    }
			}
		}
	} else {
		nlist <- make.sym.nb(nb)
		for (i in 1:n) {
			inb <- nb[[i]]
			inl <- nlist[[i]]
			if (inl[1] > 0) {
			    iwt <- wts[[i]]
			    vlist[[i]] <- numeric(length=length(inl))
			    for (j in 1:length(inl)) {
				if (inl[j] %in% inb) 
				    a <- iwt[which(inb == inl[j])]
				else a <- 0
				if (i %in% nb[[inl[j]]]) 
				    b <- wts[[inl[j]]][which(nb[[inl[j]]] == i)]
				else b <- 0
				vlist[[i]][j] <- 0.5 * (a + b)
			    }
			}
		}
	}
	res <- list(style=style, neighbours=nlist, weights=vlist)
	class(res) <- "listw"
	attr(res, "region.id") <- attr(nb, "region.id")
	attr(res, "call") <- match.call()
	attr(res, "U") <- TRUE
    }
    res
}


listw2star <- function(listw, ireg, style, n, D, a, zero.policy=NULL) {
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
    nb <- vector(mode="list", length=n)
    class(nb) <- "nb"
    wts <- vector(mode="list", length=n)
    for (i in 1:n) nb[[i]] <- 0L
    inb <- listw$neighbours[[ireg]]
    iwts <- listw$weights[[ireg]]
    cond <- TRUE
    if (inb == 0 || length(inb) == 0 || is.null(iwts)) cond <- FALSE
    if (!cond && !zero.policy) stop("No-neighbour region found")
    if (style == "W") iwts <- (n*D[ireg]*iwts) / 2
    else if (style == "S") iwts <- ((n^2)*D[ireg]*iwts) / (2*a)
    else if (style == "C") iwts <- ((n^2)*iwts) / (2*a)
    if (cond) {
    	nb[[ireg]] <- inb
    	wts[[ireg]] <- iwts
    	for (j in 1:length(inb)) {
            jj <- inb[j]
            nb[[jj]] <- ireg
            wts[[jj]] <- iwts[j]
	}
    }
    res <- list(style=style, neighbours=nb, weights=wts)
    class(res) <- c("listw", "star")
    attr(res, "region.id") <- attr(listw, "region.id")
    res
}

spdep <- function(build=FALSE) {
#	require("utils")
	.DESC <- packageDescription("spdep")
	.spdep.Version <- paste(.DESC[["Package"]], ", version ", 
		.DESC[["Version"]], ", ", .DESC[["Date"]], sep="")
	.spdep.Build <- paste("build:", .DESC[["Built"]])
	if (build) return(c(.spdep.Version, .spdep.Build))
	else return(.spdep.Version)
}
