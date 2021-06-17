# Copyright 2000-2021 by Roger S. Bivand. 
# Upgrade to sp classes February 2007
# use of dbscan 210317 #53
# s2 prototype 210612-16 not using indexing
#

dnearneigh <- function(x, d1, d2, row.names=NULL, longlat=NULL, bounds=c("GE", "LE"), use_kd_tree=TRUE, symtest=FALSE, avoid_s2=TRUE, max_cells=1000, dwithin=FALSE) {
    stopifnot(is.logical(use_kd_tree))
    use_s2_ll <- FALSE
    if (inherits(x, "SpatialPoints")) {
# correct wrong logic
        if (!is.null(longlat))
            warning("dnearneigh: longlat argument overrides object")
        if ((is.null(longlat) || !is.logical(longlat)) 
	    && !is.na(is.projected(x)) && !is.projected(x)) {
            longlat <- TRUE
        } else longlat <- FALSE
        x <- coordinates(x)
    } else {
        if (inherits(x, "sf")) {
            if (is.null(row.names)) row.names <- row.names(x)
            x <- sf::st_geometry(x)
        }
        if (inherits(x, "sfc")) {
           if (!is.null(longlat))
               warning("dnearneigh: longlat argument overrides object")
           if (!inherits(x, "sfc_POINT"))
               stop("Point geometries required")
           if (attr(x, "n_empty") > 0L) 
               stop("Empty geometries found")
           if ((is.null(longlat) || !is.logical(longlat)) 
	       && !is.na(sf::st_is_longlat(x)) && sf::st_is_longlat(x)) {
               longlat <- TRUE
           } else longlat <- FALSE
           if (longlat && sf::sf_use_s2() && !avoid_s2) {
               s2x <- sf::st_as_s2(x)
               use_s2_ll <- TRUE
           }
           x <- sf::st_coordinates(x)
        }
    }
    if (is.null(longlat) || !is.logical(longlat)) longlat <- FALSE
    if (longlat && use_kd_tree) use_kd_tree <- FALSE
    if (use_kd_tree && !requireNamespace("dbscan", quietly = TRUE)) 
        use_kd_tree <- FALSE
    if (!is.numeric(x)) stop("Data non-numeric")
    if (!is.matrix(x)) stop("Data not in matrix form")
    stopifnot(ncol(x) == 2L || ncol(x) == 3L)
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
    if (longlat && dimension > 2) stop("Only 2D spherical data accepted")
    if (!use_kd_tree && dimension > 2) stop("Only 2D without dbscan kd_tree")
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
    if (use_kd_tree) {
        z <- dbscan::frNN(x, eps=d2)$id
        z <- lapply(z, sort)
        if (d1 > 0) {
            z1 <- dbscan::frNN(x, eps=d1)$id
            z1 <- lapply(z1, sort)
            z <- lapply(seq_along(z), function(i) setdiff(z[[i]], z1[[i]])) 
        }
        z <- lapply(seq_along(z), function(i)
            {if (length(z[[i]]) == 0L) 0L else z[[i]]})
    } else if (use_s2_ll) {
        if (dwithin) {
            z <- s2::s2_dwithin_matrix(s2x, s2x, dist=d2*1000)
        } else {
            s2xb <- s2::s2_buffer_cells(s2x, distance=d2*1000,
                max_cells=max_cells)
            z <- s2::s2_intersects_matrix(s2xb, s2x)
            rm(s2xb)
        }
        z <- lapply(z, sort)
        if (d1 > 0) {
            if (dwithin) {
                z1 <- s2::s2_dwithin_matrix(s2x, s2x, dist=d1*1000)
            } else {
                s2xb <- s2::s2_buffer_cells(s2x, distance=d1*1000,
                    max_cells=max_cells)
                z1 <- s2::s2_intersects_matrix(s2xb, s2x)
                rm(s2xb)
            }
            z1 <- lapply(z1, sort)
            z <- lapply(seq_along(z), function(i) setdiff(z[[i]], z1[[i]])) 
        }
        z <- lapply(seq_along(z), function(i) setdiff(z[[i]], i))
        z <- lapply(seq_along(z), function(i)
            {if (length(z[[i]]) == 0L) 0L else z[[i]]})
    } else {
        z <- .Call("dnearneigh", d1, d2, as.integer(np), as.integer(dimension),
            x, as.integer(longlat), PACKAGE="spdep")
    }
    class(z) <- "nb"
    attr(z, "region.id") <- row.names
    attr(z, "call") <- match.call()    
    attr(z, "dnn") <- c(d1, d2)
    attr(z, "bounds") <- bounds
    attr(z, "nbtype") <- "distance"
    if (symtest) z <- sym.attr.nb(z)
    else attr(z, "sym") <- TRUE
    z
}


