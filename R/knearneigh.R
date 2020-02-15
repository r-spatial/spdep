# Copyright 2001-2014 by Roger S. Bivand. 
# Upgrade to sp classes February 2007
# Added RANN April 2010
# nn() retired 130722

knearneigh <- function(x, k=1, longlat=NULL, RANN=TRUE)
{
    if (inherits(x, "SpatialPoints")) {
        if ((is.null(longlat) || !is.logical(longlat)) 
 	   && !is.na(is.projected(x)) && !is.projected(x)) {
           longlat <- TRUE
        } else longlat <- FALSE
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
         x <- sf::st_coordinates(x)
      }
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
    if (zd) {
        if (dimension == 2) 
            warning("knearneigh: identical points found")
        else
            stop("knearneigh: identical points found")
    }
    useRANN <- RANN
    if (longlat) useRANN <- FALSE
    if (zd) useRANN <- FALSE
    if (useRANN && !requireNamespace("RANN", quietly = TRUE)) useRANN <- FALSE
# https://github.com/r-spatial/spdep/issues/38
    if (dimension > 2 && !useRANN) stop("RANN required for ncol(x) > 2")
    if (useRANN) {
#        xx <- cbind(x, out=rep(0, nrow(x)))
#        out <- as.matrix(nn(xx, p=k)$nn.idx)
# nn() retired 130722
# modified 130913 to handle zerodist points
#        zd <- zerodist(SpatialPoints(x))
#        if (nrow(zd) == 0) {
        out <- RANN::nn2(x, x, k=k+1)$nn.idx[,-1,drop=FALSE]
#        } else {
#            out0 <- nn2(x, x, k=k+1)
#            out <- t(sapply(1:np, function(i) out0$nn.idx[i,
#                -which(out0$nn.idx[i,] == i), drop=FALSE]))
#        }
        dimnames(out) <- NULL
        res <- list(nn=out, np=np, k=k, dimension=dimension, x=x)
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
    class(res) <- "knn"
    attr(res, "call") <- match.call()
    res
}
