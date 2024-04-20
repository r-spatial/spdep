# Copyright 2001-2021 by Roger S. Bivand. 
# Upgrade to sp classes February 2007
# Added RANN April 2010
# nn() retired 130722
# shift to dbscan 210317 #53
# s2 adatation 210613 https://github.com/r-spatial/s2/issues/125

knearneigh <- function(x, k=1, longlat=NULL, use_kd_tree=TRUE)
{
    use_s2_ll <- FALSE
    if (inherits(x, "SpatialPoints")) {
        if ((is.null(longlat) || !is.logical(longlat)) 
 	   && !is.na(is.projected(x)) && !is.projected(x)) {
           longlat <- TRUE
        } else longlat <- FALSE
        x <- coordinates(x)[, 1:2]
    } else if (inherits(x, "sf") || inherits(x, "sfc")) {
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
           if (longlat && sf::sf_use_s2()) {
               s2x <- sf::st_as_s2(x)
               use_s2_ll <- TRUE
           }
           x <- sf::st_coordinates(x)
        }
    } else if (inherits(x, "data.frame")) {
        x <- as.matrix(x)
    } 
    if (is.null(longlat) || !is.logical(longlat)) longlat <- FALSE
    if (!is.numeric(x)) stop("knearneigh: data non-numeric")
    if (!is.matrix(x)) stop("knearneigh: data not in matrix form")
# https://github.com/r-spatial/spdep/issues/38
    if (ncol(x) > 2 && longlat == TRUE)
        stop("ncol(x) > 2 only permitted if longlat is FALSE")
    if (any(is.na(x))) stop("knearneigh: data include NAs")
    if (longlat) {
        bb <- bbox(x)
        if (!.ll_sanity(bb))
            warning("knearneigh: coordinates are not geographical: longlat argument wrong")
    }
    if (!is.double(x)) storage.mode(x) <- "double"
    np <- nrow(x)
    dimension <- ncol(x)
# https://github.com/r-spatial/spdep/issues/38
#    if (dimension != 2) stop("knearneigh: only 2D data accepted")
    if (k >= np) stop("knearneigh: fewer data points than k")
    if (k > (np/3)) warning("k greater than one-third of the number of data points")
# modified 140117 to handle zerodist points
# (previous fix only worked for pairs and k>1)
# https://github.com/r-spatial/spdep/issues/38
    zd <- any(duplicated(x)) #zerodist(SpatialPoints(x))
# kNN() handles duplicate points
    if (zd && !use_s2_ll)  {
        warning("knearneigh: identical points found")
        if (use_kd_tree) {
            warning("knearneigh: kd_tree not available for identical points")
            use_kd_tree <- FALSE
        }
        if (dimension > 2) stop("kd_tree required for more than 2 dimensions")
    }
    if (longlat) use_kd_tree <- FALSE
    if (use_kd_tree && !requireNamespace("dbscan", quietly = TRUE)) use_kd_tree <- FALSE
# https://github.com/r-spatial/spdep/issues/38
    
    if (dimension > 2 && use_s2_ll) use_s2_ll <- FALSE
    if (use_kd_tree) {
        out <- dbscan::kNN(x, k=k)$id
        dimnames(out) <- NULL
        res <- list(nn=out, np=np, k=k, dimension=dimension, x=x)
    } else {
        if (use_s2_ll) {
            z <- s2::s2_closest_edges(s2x, s2x, k=(k+1))
            z <- lapply(z, sort)
            z <- lapply(seq_along(z), function(i) setdiff(z[[i]], i))
            res <- list(nn=do.call("rbind", z), np=np, k=k,
    	        dimension=dimension, x=x)
        } else {
            xx <- c(x[,1], x[,2])
            storage.mode(xx) <- "double"
            nn <- integer(np*k)
            dnn <- double(np*k)
            z <- .C("knearneigh", k=as.integer(k), np=as.integer(np),
                dimension=as.integer(dimension),
                xx=xx, nn=as.integer(nn), dnn=dnn,
	        as.integer(longlat), PACKAGE="spdep")
            res <- list(nn=matrix(z$nn, np, k, byrow=TRUE), np=np, k=k,
    	        dimension=dimension, x=x)
        }
    }
    class(res) <- "knn"
    attr(res, "call") <- match.call()
    res
}
