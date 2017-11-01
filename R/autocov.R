# Calculates the autocovariate to be used in autonormal, autopoisson 
# or autologistic regression. Three distance-weighting schemes are available
#
# z is the response variable
# xy is the matrix of coordinates
# nbs is "neighbourhood size", selected by user; default is 1
# type defines the weighting scheme: 
#	"one" gives equal weight to all data points in the neighbourhood; 
#	"inverse" (the default) weights by inverse Euclidean distance;
#	"inverse.squared" weights by the square of "inverse"
#
# by Carsten F. Dormann, 04.11.2005, carsten.dormann@ufz.de
# Re-implementation allowing list representation
# Roger Bivand 28.11.2005
# Upgrade to sp classes February 2007, longlat sanity June 2010
# style changed from "W" to "B" 2015-01-27; see Bardos et al. (2015) for details

autocov_dist <- function(z, xy, nbs=1, type="inverse", zero.policy=NULL,
   style="B", longlat=NULL) {
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
    stopifnot(is.vector(z))
   if (type=="one") expo <- 0
   if (type=="inverse") expo <- 1
   if (type=="inverse.squared") expo <- 2
   if (inherits(xy, "SpatialPoints")) {
      if ((is.null(longlat) || !is.logical(longlat)) 
	 && !is.na(is.projected(xy)) && !is.projected(xy)) {
         longlat <- TRUE
      } else longlat <- FALSE
      xy <- coordinates(xy)
   } else if (is.null(longlat) || !is.logical(longlat)) longlat <- FALSE
   stopifnot(ncol(xy) == 2)
   if (longlat) {
        bb <- bbox(xy)
        if (!.ll_sanity(bb))
            warning("Coordinates are not geographical: longlat argument wrong")
   }
   nb <- dnearneigh(xy, 0, nbs, longlat=longlat)
   if (any(card(nb) == 0)) warning(paste("With value", nbs,
      "some points have no neighbours"))
   nbd <- nbdists(nb, xy, longlat=longlat)
   if (expo == 0) lw <- nb2listw(nb, style=style, zero.policy=zero.policy)
   else {
      gl <- lapply(nbd, function(x) 1/(x^expo))
      lw <- nb2listw(nb, glist=gl, style=style, zero.policy=zero.policy)
   }
   lag(lw, z, zero.policy=zero.policy)
}

.ll_sanity <- function(bb) {
        TOL <- get_ll_TOL()
	tol <- .Machine$double.eps ^ TOL
	W <- bb[1,1] < -180 && 
	    !isTRUE(all.equal((bb[1, 1] - -180), 0, tolerance = tol))
        if (W) attr(W, "out") <- bb[1,1]
	E <- bb[1,2] > 360 && 
	    !isTRUE(all.equal((bb[1, 2] - 360), 0, tolerance = tol))
        if (E) attr(E, "out") <- bb[1,2]
	S<- bb[2,1] < -90 && 
	    !isTRUE(all.equal((bb[2, 1] - -90), 0, tolerance = tol))
        if (S) attr(S, "out") <- bb[2,1]
	N <- bb[2,2] > 90 && 
	    !isTRUE(all.equal((bb[2, 2] - 90), 0, tolerance = tol))
        if (N) attr(N, "out") <- bb[2,2]
        res <- !(any(W || E || S || N))
        attr(res, "details") <- list(W, E, S, N)
	res
}

