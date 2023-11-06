# Copyright 2001-10 by Roger S. Bivand and Virgilio Gomez-Rubio
#

nb2listw <- function(neighbours, glist=NULL, style="W", zero.policy=NULL)
{
	if(!inherits(neighbours, "nb")) stop("Not a neighbours list")
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
	if (!(style %in% c("W", "B", "C", "S", "U", "minmax")))
		stop(paste("Style", style, "invalid"))
	n <- length(neighbours)
	if (n < 1) stop("non-positive number of entities")
	cardnb <- card(neighbours)
	if (!zero.policy)
		if (any(cardnb == 0)) stop("Empty neighbour sets found")
	vlist <- vector(mode="list", length=n)
	if (is.null(glist)) {
		glist <- vector(mode="list", length=n)
		for (i in 1:n)
			if(cardnb[i] > 0) {
				glist[[i]] <- rep(1, cardnb[i])
				mode(glist[[i]]) <- "numeric"
			}
		attr(vlist, "mode") <- "binary"
	} else {
		if (length(glist) != n) stop("glist wrong length")
		if (any(cardnb != unlist(lapply(glist, length))))
			stop("neighbours and glist do not conform")
		if (any(is.na(unlist(glist))))
			stop ("NAs in general weights list")
		if (!zero.policy && any(sapply(glist, function(x) 
			isTRUE(all.equal(sum(x), 0)))))
			warning("zero sum general weights") 
		glist <- lapply(glist, function(x) {mode(x) <- "numeric"; x})
		attr(vlist, "mode") <- "general"
		attr(vlist, "glist") <- glist
		attr(vlist, "glistsym") <- is.symmetric.glist(neighbours, glist)
	}
	attr(vlist, as.character(style)) <- TRUE
	if (zero.policy) {
		eff.n <- n - length(which(cardnb == 0))
		if (eff.n < 1) stop("No valid observations")
	} else eff.n <- n
	if (style == "W") {
		d <- unlist(lapply(glist, sum))
		for (i in 1:n) {
			if (cardnb[i] > 0) {
			    if (d[i] > 0) vlist[[i]] <- (1/d[i]) * glist[[i]]
			    else vlist[[i]] <- 0 * glist[[i]]
			}
		}
		attr(vlist, "comp") <- list(d=d)
	}
	if (style == "B") {
		for (i in 1:n) {
			if (cardnb[i] > 0) vlist[[i]] <- glist[[i]]
		}
	}
	if (style == "C" || style == "U" || style == "minmax") {
		D <- sum(unlist(glist))
		if (is.na(D) || !(D > 0))
			stop(paste("Failure in sum of weights:", D))
		for (i in 1:n) {
			if (cardnb[i] > 0) {
				if (style == "C")
					vlist[[i]] <- (eff.n/D) * glist[[i]]
				else if(style == "U")
                                        vlist[[i]] <- (1/D) * glist[[i]]
                                else vlist[[i]] <- glist[[i]]
			}
		}
	}
	if (style == "S") {
		glist2 <- lapply(glist, function(x) x^2)
		q <- sqrt(unlist(lapply(glist2, sum)))
		for (i in 1:n) {
			if (cardnb[i] > 0) {
			    if (q[i] > 0) glist[[i]] <- (1/q[i]) * glist[[i]]
			    else glist[[i]] <- 0 * glist[[i]]
			}
		}
		Q <- sum(unlist(glist))
		if (is.na(Q) || !(Q > 0))
		    stop(paste("Failure in sum of intermediate weights:", Q))
		for (i in 1:n) {
			if (cardnb[i] > 0)
				vlist[[i]] <- (eff.n/Q) * glist[[i]]
		}
		attr(vlist, "comp") <- list(q=q, Q=Q, eff.n=eff.n)
	}
	style <- style
	if (!zero.policy)
		if (any(is.na(unlist(vlist))))
			stop ("NAs in coding scheme weights list")
        if (style == "minmax") {
            res <- list(style=style, neighbours=neighbours, weights=vlist)
	    class(res) <- c("listw", "nb")
            mm <- minmax.listw(res)
            vlist <- lapply(vlist, function(x) (1/c(mm)) * x)
        }
	res <- list(style=style, neighbours=neighbours, weights=vlist)
	class(res) <- c("listw", "nb")
	attr(res, "region.id") <- attr(neighbours, "region.id")
	attr(res, "call") <- match.call()
        attr(res, "zero.policy") <- zero.policy
	if (!is.null(attr(neighbours, "GeoDa")))
		attr(res, "GeoDa") <- attr(neighbours, "GeoDa")
	if (!is.null(attr(res, "GeoDa")$dist)) 
		attr(res, "GeoDa")$dist <- NULL
	res
}


nb2WB <- function(nb)
{
# class to inherits Jari Oksanen 080603
  	if (!inherits(nb, "nb")) stop("not a neighbours list")
        num <- card(nb)
        if (any(num == 0)) nb[num == 0] <- NULL
        adj <- unlist(nb)
        weights <- rep(1, sum(num))

        list(adj=adj, weights=weights, num=num)
}

listw2WB <- function(listw)
{
	if (!inherits(listw, "listw")) stop("not listw class object")
        num <- card(listw$neighbours)
        if (any(num == 0)) listw$neighbours[num == 0] <- NULL
        adj <- unlist(listw$neighbours)
        weights <- unlist(listw$weights)

        list(adj=adj, weights=weights, num=num)
}

minmax.listw <- function(listw) {
    W <- listw2mat(listw)
    rm <- max(rowSums(W))
    cm <- max(colSums(W))
    res <- min(c(rm, cm))
    attr(res, "rowmax") <- rm
    attr(res, "colmax") <- cm
    res
}

