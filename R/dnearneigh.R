# Copyright 2000-2014 by Roger S. Bivand. 
# Upgrade to sp classes February 2007
#

dnearneigh <- function(x, d1, d2, row.names=NULL, longlat=NULL, bounds=c("GT", "LE")) {
   if (inherits(x, "SpatialPoints")) {
# correct logic
      if (!is.null(longlat))
          warning("dnearneigh: longlat overriden for Spatial object")
      if (!is.na(is.projected(x)) && !is.projected(x)) {
          longlat <- TRUE
      } else longlat <- FALSE
      x <- coordinates(x)[, 1:2]
   } else if (is.null(longlat) || !is.logical(longlat)) longlat <- FALSE
    if (!is.numeric(x)) stop("Data non-numeric")
    if (!is.matrix(x)) stop("Data not in matrix form")
    stopifnot(ncol(x) == 2)
    if (any(is.na(x))) stop("Data include NAs")
    if (longlat) {
        bb <- bbox(x)
        if (!.ll_sanity(bb))
            warning("Coordinates are not geographical: longlat argument wrong")
    }
#    if (!is.double(x)) storage.mode(x) <- "double"
    np <- nrow(x)
    if (np < 1) stop("non-positive number of rows in x")
    if (!is.null(row.names)) {
	if(length(row.names) != np)
            stop("row.names wrong length")
	if (length(unique(row.names)) != length(row.names))
	    stop("non-unique row.names given")
    }
    if (is.null(row.names)) row.names <- as.character(1:np)
    dimension <- ncol(x)
    if (dimension > 2) stop("Only 2D data accepted")
    md <- 0
    if (d1 < 0) d1 <- 0.0
    if (!longlat) {
	for (i in 1:dimension) md <- sum(md, (diff(range(x[,i]))^2))
	md <- md + (.Machine$double.eps)^(1/4)
    	if (d2 > sqrt(md)) d2 <- sqrt(md)
    }
    stopifnot(is.character(bounds))
    stopifnot(length(bounds) == 2)
    stopifnot(isTRUE(bounds[1] %in% c("GE", "GT")))
    stopifnot(isTRUE(bounds[2] %in% c("LE", "LT")))
    storage.mode(x) <- "double"
    storage.mode(d1) <- "double"
    storage.mode(d2) <- "double"
    attr(d1, "equal") <- bounds[1] == "GE"
    attr(d2, "equal") <- bounds[2] == "LE"
    z <- .Call("dnearneigh", d1, d2, as.integer(np),
        as.integer(dimension), x, as.integer(longlat), 
	PACKAGE="spdep")
    attr(z[[1]], "region.id") <- row.names
    attr(z[[1]], "call") <- match.call()
    attr(z[[1]], "dnn") <- c(d1, d2)
    attr(z[[1]], "bounds") <- bounds
    z[[1]] <- sym.attr.nb(z[[1]])
    z[[1]]
}
