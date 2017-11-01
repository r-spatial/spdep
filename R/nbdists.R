# Copyright 2001-2010 by Roger Bivand
# Upgrade to sp classes February 2007
#


nbdists <- function(nb, coords, longlat=NULL) {
	if (!inherits(nb, "nb")) 
        	stop("Not a neighbours list")
   	if (inherits(coords, "SpatialPoints")) {
      		if ((is.null(longlat) || !is.logical(longlat)) 
		    && !is.na(is.projected(coords)) && !is.projected(coords)) {
         		longlat <- TRUE
      		} else longlat <- FALSE
      		coords <- coordinates(coords)[, 1:2]
   	} else if (is.null(longlat) || !is.logical(longlat)) longlat <- FALSE
	if (!is.numeric(coords)) stop("Data non-numeric")
	if (!is.matrix(coords)) 
            stop("Data not in matrix form")
        stopifnot(ncol(coords) == 2)
        if (any(is.na(coords))) 
            stop("Data include NAs")
        if (longlat) {
            bb <- bbox(coords)
            if (!.ll_sanity(bb))
                warning("Coordinates are not geographical: longlat argument wrong")
        }
	if (!is.double(coords)) storage.mode(coords) <- "double"
	n.nb <- length(nb)
	np <- nrow(coords)
        if (np != n.nb) 
            stop("Number of coords not equal to number of regions")
        dimension <- ncol(coords)
        dlist <- .Call("nbdists", nb, as.matrix(coords), as.integer(np), 
            as.integer(dimension), as.integer(longlat), PACKAGE="spdep")
	attr(dlist[[1]], "call") <- match.call()
	dlist[[1]]
}

