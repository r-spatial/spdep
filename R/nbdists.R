# Copyright 2001-2019 by Roger Bivand
# Upgrade to sp classes February 2007
# s2 prototype 210612
#


nbdists <- function(nb, coords, longlat=NULL) {
	if (!inherits(nb, "nb")) 
        	stop("Not a neighbours list")
        use_s2_ll <- FALSE
   	if (inherits(coords, "SpatialPoints")) {
                if (!is.null(longlat))
                    warning("dnearneigh: longlat argument overrides object")
      		if ((is.null(longlat) || !is.logical(longlat)) 
		    && !is.na(is.projected(coords)) && !is.projected(coords)) {
         		longlat <- TRUE
      		} else longlat <- FALSE
      		coords <- coordinates(coords)[, 1:2]
        } else {
            if (inherits(coords, "sf")) {
                if (is.null(row.names)) row.names <- row.names(coords)
                coords <- sf::st_geometry(coords)
            }
            if (inherits(coords, "sfc")) {
                if (!is.null(longlat))
                    warning("dnearneigh: longlat argument overrides object")
                if (!inherits(coords, "sfc_POINT"))
                    stop("Point geometries required")
                if (attr(coords, "n_empty") > 0L) 
                    stop("Empty geometries found")
                if ((is.null(longlat) || !is.logical(longlat)) 
	            && !is.na(sf::st_is_longlat(coords)) && 
                    sf::st_is_longlat(coords)) {
                    longlat <- TRUE
                } else longlat <- FALSE
                if (longlat && sf::sf_use_s2()) {
                    s2x <- sf::st_as_s2(coords)
                    use_s2_ll <- TRUE
                }
                coords <- sf::st_coordinates(coords)
            }
        }
        if (is.null(longlat) || !is.logical(longlat)) longlat <- FALSE
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
        if (use_s2_ll) {
            dlist <- vector(mode="list", length=1L)
            dlist[[1]] <- lapply(seq_along(nb), function(i) {
                d_i <- s2::s2_distance(s2x[i], s2x[nb[[i]]])
                units::set_units(units::set_units(d_i, "m"), "km")})
        } else {
            dlist <- .Call("nbdists", nb, as.matrix(coords), as.integer(np), 
                as.integer(dimension), as.integer(longlat), PACKAGE="spdep")
        }
	attr(dlist[[1]], "call") <- match.call()
	dlist[[1]]
}

