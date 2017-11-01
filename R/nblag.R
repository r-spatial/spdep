# Copyright 2001-2006 by Roger Bivand
#


nblag <- function(neighbours, maxlag)
{
	if (!is.null(attributes(neighbours)$self.included) &&
		(as.logical(attributes(neighbours)$self.included)))
		stop("No lags for neighbours lists including self")
	n <- length(neighbours)
	if (n < 1) stop("non-positive number of entities")
	if (maxlag < 2) stop("maxlag less than 2")
	lags <- vector(mode="list", length=maxlag)
	lags[[1]] <- neighbours
	cds <- card(neighbours)
	for (thislag in 2:maxlag)
		lags[[thislag]] <- vector(mode="list", length=n)
	for (i in 1:n) {
		already <- i
		new <- neighbours[[i]]
		for (thislag in 2:maxlag) {
			if (cds[i] > 0) {
				already <- c(already, new)
				active <- new
				new <- NULL
				for (j in active)
					new <- c(new, neighbours[[j]])
				new <- sort(unique(new))
				res <- new[-which(new %in% already)]
				if(length(res) == 0L) 
					lags[[thislag]][[i]] <- 0L
				else lags[[thislag]][[i]] <- res
			}
			else lags[[thislag]][[i]] <- 0L
		}
	}
	for (i in 1:maxlag) {
		class(lags[[i]]) <- "nb"
		attr(lags[[i]], "region.id") <- attr(neighbours, "region.id")
		lags[[i]] <- sym.attr.nb(lags[[i]])
	}
	attr(lags, "call") <- match.call()
	lags
}

# Copyright 2006-2010 (c) Giovanni Millo and Roger Bivand

nblag_cumul <- function (nblags) {
    if (any(sapply(nblags, function(x) class(x) != "nb")))
        stop("nblags must be a list of neighbour objects")
    maxlag <- length(nblags)
    if (maxlag < 2) stop("maxlag must be greater than 1")

    n <- length(nblags[[1]])
    lags <- vector(mode="list", length=n)
    ncard <- card(nblags[[1]])
    neigh <- which(ncard > 0)
    nneigh <- which(ncard == 0)
    for (i in nneigh) lags[[i]] <- 0L
    for (i in neigh) {
        res <- nblags[[1]][[i]]
	for (j in 2:maxlag) res <- c(res, nblags[[j]][[i]])
        res <- as.integer(res[order(unique(res))])
        if (any(res == 0)) res <- res[-which(res == 0)]
        lags[[i]] <- res
    }
    attr(lags, "region.id") <- attr(nblags[[1]], "region.id")
    attr(lags, "call") <- match.call()
    class(lags) <- "nb"
    lags
}
