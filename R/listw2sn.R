# Copyright 2001-7 by Roger Bivand
#

as.spam.listw <- function(listw) {
    if (requireNamespace("spam", quietly = TRUE)) {
#if (!require(spam)) stop("spam not available")
        N <- length(listw$neighbours)
        W_sn <- listw2sn(listw)
        rpts <- as.integer(cumsum(c(1, card(listw$neighbours))))
        W <- new("spam", entries=W_sn$weights, colindices=W_sn$to,
            rowpointers=rpts, dimension=as.integer(c(N, N)))
        stopifnot(spam::validspamobject(W))
        return(W)
    } else stop("spam not available")
}

listw2U_spam <- function(lw) 0.5 * (lw + t(lw))


listw2sn <- function(listw) {
	if(!inherits(listw, "listw")) stop("not a listw object")
	cardw <- card(listw$neighbours)
	scard <- sum(cardw)
	z <- .Call("listw2sn", listw$neighbours, listw$weights,
		as.integer(cardw), as.integer(scard), PACKAGE="spdep")
	res <- as.data.frame(list(from=z[[1]], to=z[[2]], weights=z[[3]]))
	class(res) <- c(class(res), "spatial.neighbour")
	attr(res, "n") <- length(attr(listw, "region.id"))
	attr(res, "region.id") <- attr(listw, "region.id")
	neighbours.attrs <- names(attributes(listw$neighbours))
	attr(res, "neighbours.attrs") <- neighbours.attrs
	weights.attrs <- names(attributes(listw$weights))
	attr(res, "weights.attrs") <- weights.attrs
	if (!(is.null(attr(listw, "GeoDa"))))
		attr(res, "GeoDa") <- attr(listw, "GeoDa")
	attr(res, "listw.call") <- attr(listw, "call")
	res
}

sn2listw <- function(sn) {
	if(!inherits(sn, "spatial.neighbour")) 
	    stop("not a spatial.neighbour object")
	n <- attr(sn, "n")
	if (n < 1) stop("non-positive n")
	region.id <- attr(sn, "region.id")
        stopifnot(all(!is.na(sn[,1])))
        stopifnot(all(!is.na(sn[,2])))
        stopifnot(all(!is.na(sn[,3])))
	nlist <- vector(mode="list", length=n)
	class(nlist) <- "nb"
	attr(nlist, "region.id") <- region.id
	vlist <- vector(mode="list", length=n)
	rle.sn <- rle(sn[,1])
	if (n != length(rle.sn$lengths)) {
            nnhits <- region.id[which(!(1:n %in% rle.sn$values))]
	    warning(paste(paste(nnhits, collapse=", "),
                ifelse(length(nnhits) < 2, "is not an origin",
                "are not origins")))
        }
	cs1.sn <- cumsum(rle.sn$lengths)
	cs0.sn <- c(1, cs1.sn[1:(n-1)]+1)
	ii <- 1
	for (i in 1:n) {
		if (!is.na(rle.sn$value[ii]) && rle.sn$value[ii] == i) {
                        ni <- as.integer(sn[cs0.sn[ii]:cs1.sn[ii],2])
                        o <- order(ni)
			nlist[[i]] <- ni[o]
			vlist[[i]] <- as.double(sn[cs0.sn[ii]:cs1.sn[ii],3])[o]
			ii <- ii+1
		} else {
			nlist[[i]] <- 0L
		}
	}
	res <- list(style=as.character(NA), neighbours=nlist, weights=vlist)
	class(res) <- c("listw", "nb")
	if (!(is.null(attr(sn, "GeoDa"))))
		attr(res, "GeoDa") <- attr(sn, "GeoDa")
	attr(res, "region.id") <- region.id
	attr(res, "call") <- match.call()
	res
}

