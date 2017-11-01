# Copyright 2004-2013 by Roger Bivand 
#

nb2blocknb <- function(nb=NULL, ID, row.names = NULL) {
        # Jacquelyn Pless suggestion 131204
        if (is.null(nb)) {
            blks <- unique(as.character(ID))
            nb <- lapply(blks, function(x) 0L)
            class(nb) <- "nb"
            attr(nb, "region.id") <- blks
        }
	if (!inherits(nb, "nb")) stop("not an nb object")
	nbNames <- as.character(attr(nb, "region.id"))
	entNames <- as.character(ID)
	if (!identical(sort(nbNames), sort(unique(entNames))))
		stop("names do not match exactly")
	n <- length(entNames)
	if (n < 1) stop("non-positive number of entities")
	if (!is.null(row.names)) {
		if (length(row.names) != n) 
			stop("row.names wrong length")
		if (length(unique(row.names)) != length(row.names)) 
		stop("non-unique row.names given")
	} else {
		row.names <- as.character(1:n)
	}
	inter <- lapply(as.list(nbNames), 
		function(x) which(match(entNames, x) == 1))

	res <- vector(mode="list", length=n)
	for (i in 1:n) {
		ii <- match(entNames[i], nbNames)
		blocks <- c(ii, nb[[ii]])
		vec <- sort(unlist(inter[blocks]))
                svec <- vec[vec != i]
                if (length(svec) == 0) {
                    res[[i]] <- 0L
                } else {
                    res[[i]] <- svec
                }
# Ann Hartell bug for NULL 2016-08-22
#		res[[i]] <- ifelse(length(svec) == 0, 0L, svec)
	}

	attr(res, "region.id") <- row.names
	class(res) <- "nb"
	attr(res, "block") <- TRUE
	attr(res, "call") <- match.call()
	res <- sym.attr.nb(res)
	res
}


