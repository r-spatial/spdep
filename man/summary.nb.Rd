% Copyright 2001-19 by Roger S. Bivand
\name{summary.nb}
\alias{summary.nb}
\alias{print.summary.nb}
\alias{print.nb}
\alias{summary.listw}
\alias{print.summary.listw}
\alias{print.listw}
\title{Print and summary function for neighbours and weights lists}
\description{
The function prints summary measures for links in a neighbours list. If a
matrix of coordinates is given as well, summary descriptive measures for
the link lengths are also printed. Print and summary functions are also available for \code{"listw"} weights list objects, also reporting constants (S0, S1, S2) used in inference for global spatial autocorrelation statistics such as Moran's I, Geary's C, join-count tests and Getis-Ord G.
}
\usage{
\method{summary}{nb}(object, coords=NULL, longlat = NULL, scale = 1, ...)
\method{print}{nb}(x, ...)
\method{summary}{listw}(object, coords, longlat, zero.policy = attr(object, "zero.policy"), 
 scale = 1, adjust.n=TRUE, ...)
\method{print}{listw}(x, zero.policy = attr(x, "zero.policy"), ...)
}

\arguments{
  \item{object}{an object of class \code{nb}}
  \item{coords}{matrix of region point coordinates or a SpatialPoints object or an \code{sfc} points object}
  \item{longlat}{TRUE if point coordinates are longitude-latitude decimal degrees, in which case distances are measured in kilometers; if coords is a SpatialPoints object, the value is taken from the object itself}
  \item{...}{additional arguments affecting the output produced}
  \item{x}{an object of class \code{nb}}
  \item{zero.policy}{default \code{attr(listw, "zero.policy")} as set when \code{listw} was created, if attribute not set, use global option value; if FALSE stop with error for any empty neighbour sets}
  \item{scale}{passed through to \code{stem()} for control of plot length}
  \item{adjust.n}{default TRUE, if FALSE the number of observations is not adjusted for no-neighbour observations, if TRUE, the number of observations in \code{spweights.constants} is adjusted}
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{plot.nb}}}

\examples{
columbus <- st_read(system.file("shapes/columbus.shp", package="spData")[1], quiet=TRUE)
col.gal.nb <- read.gal(system.file("weights/columbus.gal", package="spData")[1])
coords <- st_centroid(st_geometry(columbus), of_largest_polygon=TRUE)
col.gal.nb
summary(col.gal.nb, coords)
col.listw <- nb2listw(col.gal.nb, style="W")
col.listw
summary(col.listw)
col_geoms <- st_geometry(columbus)
col_geoms[21] <- st_buffer(col_geoms[21], dist=-0.05)
st_geometry(columbus) <- col_geoms
nb <- poly2nb(columbus)
summary(nb)
try(nb2listw(nb, style="W"))
summary(nb2listw(nb, style="W", zero.policy=TRUE))
summary(nb2listw(nb, style="W", zero.policy=TRUE), adjust.n=FALSE)
}
\keyword{spatial}
