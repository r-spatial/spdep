# Copyright 2001-10 by Roger Bivand, Markus Reder and Werner Mueller, 2015 Martin Gubri
#


nb2mat <- function(neighbours, glist=NULL, style="W", zero.policy=NULL)
{
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
	if(!inherits(neighbours, "nb")) stop("Not a neighbours list")
	listw <- nb2listw(neighbours, glist=glist, style=style,
		zero.policy=zero.policy)
	res <- listw2mat(listw)
	attr(res, "call") <- match.call()
	res
}

listw2mat <- function(listw) {
	n <- length(listw$neighbours)
	if (n < 1) stop("non-positive number of entities")
	cardnb <- card(listw$neighbours)
	if (any(is.na(unlist(listw$weights))))
		stop ("NAs in general weights list")
	res <- matrix(0, nrow=n, ncol=n)
	for (i in 1:n)
	    if (cardnb[i] > 0)
		res[i, listw$neighbours[[i]]] <- listw$weights[[i]]
	if (!is.null(attr(listw, "region.id")))
		row.names(res) <- attr(listw, "region.id")
	res
}


mat2listw <- function(x, row.names=NULL, style="M") {
	if (!(is.matrix(x) || is(x, "sparseMatrix"))) stop("x is not a matrix")
	n <- nrow(x)
	if (n < 1) stop("non-positive number of entities")
	m <- ncol(x)
	if (n != m) stop("x must be a square matrix")
	if (any(x < 0)) stop("values in x cannot be negative")
	if (any(is.na(x))) stop("NA values in x not allowed")
    	if (!is.null(row.names)) {
		if(length(row.names) != n)
            		stop("row.names wrong length")
		if (length(unique(row.names)) != length(row.names))
	    		stop("non-unique row.names given")
    	}
    	if (is.null(row.names)) {
		if (!is.null(row.names(x))) {
			row.names <- row.names(x)
		} else {
			row.names <- as.character(1:n)
		}
	}
#	style <- "M"
        if (is(x, "sparseMatrix")) {
            xC <- as(x, "dgCMatrix")
            i <- slot(xC, "i")+1
            p <- slot(xC, "p")
            dp <- diff(p)
            rp <- rep(seq_along(dp), dp)
            df0 <- data.frame(from=i, to=rp, weights=slot(xC, "x"))
            o <- order(df0$from, df0$to)
            df <- df0[o,]
            class(df) <- c(class(df), "spatial.neighbour")
            attr(df, "region.id") <- row.names
            attr(df, "n") <- dim(xC)[1]
            res <- sn2listw(df)
            neighbours <- res$neighbours
            weights <- res$weights
        } else {
	    neighbours <- vector(mode="list", length=n)
	    weights <- vector(mode="list", length=n)
	    for (i in 1:n) {
		nbs  <- which(x[i,] > 0.0)
		if (length(nbs) > 0) {
			neighbours[[i]] <- nbs
			weights[[i]] <- as.double(x[i, nbs]) # Laurajean Lewis
		} else {
			neighbours[[i]] <- 0L
		}
	    }
        }
	attr(weights, "mode") <- "unknown" # Brian Rubineau
	class(neighbours) <- "nb"
	attr(neighbours, "region.id") <- row.names
 	attr(neighbours, "call") <- NA
        attr(neighbours, "sym") <- is.symmetric.nb(neighbours, 
		verbose=FALSE, force=TRUE)
	res <- list(style=style, neighbours=neighbours, weights=weights)
	class(res) <- c("listw", "nb")
	attr(res, "region.id") <- attr(neighbours, "region.id")
	attr(res, "call") <- match.call()
        if (style != "M") {
            res <- nb2listw(res$neighbours, glist=res$weights, style=style,
                zero.policy=TRUE)
        }
	res
}
