# Copyright 2000-2021 by Roger S. Bivand. 
# Upgrade to sp classes February 2007
#

dnearneigh <- function(x, d1, d2, row.names=NULL, longlat=NULL, bounds=c("GE", "LE"), prefer_heuristic=TRUE, heuristic_n=1000, heuristic_fuzz=0.02, rtree=FALSE, exact_heuristic=FALSE, symtest=FALSE) {
    heur_OK <- FALSE
# rtree code OP in rtree branch
#    if (rtree) rtree <- FALSE        
    if (inherits(x, "SpatialPoints")) {
# correct wrong logic
        if (!is.null(longlat))
            warning("dnearneigh: longlat argument overrides object")
        if ((is.null(longlat) || !is.logical(longlat)) 
	    && !is.na(is.projected(x)) && !is.projected(x)) {
            longlat <- TRUE
        } else longlat <- FALSE
        xc <- st_as_sfc(as(x, "SpatialPoints"))
        heur_OK <- TRUE
        x <- coordinates(x)[, 1:2]
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
           xc <- x
           heur_OK <- TRUE
           x <- sf::st_coordinates(xc)
        }
    }
    if (is.null(longlat) || !is.logical(longlat)) longlat <- FALSE
    if (longlat) heur_OK <- FALSE
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
    stopifnot(is.logical(prefer_heuristic))
    if (np > heuristic_n && prefer_heuristic && heur_OK) {
      if (exact_heuristic) {
        if (rtree) {
            xtree <- rtree::RTree(x)
            outer_d2_buf_ints <- rtree::withinDistance.RTree(xtree, x,
                d2+heuristic_fuzz)
        } else {
            outer_d2_buf <- sf::st_buffer(xc, d2*(1+heuristic_fuzz))
            outer_d2_buf_ints <- sf::st_intersects(outer_d2_buf, xc)
        }
        outer_d2_buf_ints1 <- lapply(seq_along(outer_d2_buf_ints),
            function(i) {outer_d2_buf_ints[[i]][outer_d2_buf_ints[[i]] != i]})
        if (rtree) {
            inner_d2_buf_ints <- rtree::withinDistance.RTree(xtree, x, 
                d2-heuristic_fuzz)
        } else {
            inner_d2_buf <- sf::st_buffer(xc, d2*(1-heuristic_fuzz))
            inner_d2_buf_ints <- sf::st_intersects(inner_d2_buf, xc)
        }
        inner_d2_buf_ints1 <- lapply(seq_along(inner_d2_buf_ints),
            function(i) {inner_d2_buf_ints[[i]][inner_d2_buf_ints[[i]] != i]})
        inner_d2_buf_diff <- lapply(seq_along(outer_d2_buf_ints1), 
            function(i) setdiff(outer_d2_buf_ints1[[i]],
                inner_d2_buf_ints1[[i]]))
        z2 <- .Call("dnearneigh1", d1, d2, as.integer(np), x,
            inner_d2_buf_diff, PACKAGE="spdep")
        z <- lapply(seq_along(z2),
            function(i) sort(union(inner_d2_buf_ints1[[i]], z2[[i]])))
        if (d1 > 0.0) {
            if (rtree) {
                outer_d1_buf_ints <- rtree::withinDistance.RTree(xtree, x, 
                    d1+heuristic_fuzz)
            } else {
                outer_d1_buf <- sf::st_buffer(xc, d1*(1+heuristic_fuzz))
                outer_d1_buf_ints <- sf::st_intersects(outer_d1_buf, xc)
            }
            outer_d1_buf_ints1 <- lapply(seq_along(outer_d1_buf_ints),
                function(i) {outer_d1_buf_ints[[i]][outer_d1_buf_ints[[i]] != i]})
            d2_outer_d1 <- lapply(seq_along(z), 
                function(i) setdiff(z[[i]], outer_d1_buf_ints1[[i]]))
            if (rtree) {
                inner_d1_buf_ints <- rtree::withinDistance.RTree(xtree, x, 
                    d1-heuristic_fuzz)
            } else {
                inner_d1_buf <- sf::st_buffer(xc, d1*(1-heuristic_fuzz))
                inner_d1_buf_ints <- sf::st_intersects(inner_d1_buf, xc)
            }
            inner_d1_buf_ints1 <- lapply(seq_along(inner_d1_buf_ints),
                function(i) {inner_d1_buf_ints[[i]][inner_d1_buf_ints[[i]] != i]})
            inner_d1_buf_diff <- lapply(seq_along(outer_d1_buf_ints1), 
                function(i) setdiff(outer_d1_buf_ints1[[i]],
                    inner_d1_buf_ints1[[i]]))
            z1 <- .Call("dnearneigh1", d1, d2, as.integer(np), x,
                inner_d1_buf_diff, PACKAGE="spdep")
            z <- lapply(seq_along(z1),
                function(i) sort(union(d2_outer_d1[[i]], z1[[i]])))
        }
      } else {
        if (rtree) {
            xtree <- rtree::RTree(x)
            z <- rtree::withinDistance.RTree(xtree, x, d2)
            z <- lapply(seq_along(z), function(i) {z[[i]][z[[i]] != i]})
            if (d1 > 0) {
                z1 <- rtree::withinDistance.RTree(xtree, x, d1)
                z1 <- lapply(seq_along(z1), function(i) {z1[[i]][z1[[i]] != i]})
                z <- lapply(seq_along(z), function(i) setdiff(z[[i]], z1[[i]])) 
            }
        } else {
            z_buf <- sf::st_buffer(xc, d2)
            z <- sf::st_intersects(z_buf, xc)
            z <- lapply(seq_along(z), function(i) {z[[i]][z[[i]] != i]})
            if (d1 > 0) {
                z1_buf <- sf::st_buffer(xc, d1)
                z1 <- sf::st_intersects(z1_buf, xc)
                z1 <- lapply(seq_along(z1), function(i) {z1[[i]][z1[[i]] != i]})
                z <- lapply(seq_along(z), function(i) setdiff(z[[i]], z1[[i]])) 
            }
        }
      }
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


