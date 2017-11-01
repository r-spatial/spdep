% Copyright 2001-3 by Roger S. Bivand
\name{summary.nb}
\alias{summary.nb}
\alias{print.nb}
\alias{summary.listw}
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
\method{summary}{listw}(object, coords, longlat, zero.policy = NULL,
 scale = 1, ...)
\method{print}{listw}(x, zero.policy = NULL, ...)
}

\arguments{
  \item{object}{an object of class \code{nb}}
  \item{coords}{matrix of region point coordinates or a SpatialPoints object}
  \item{longlat}{TRUE if point coordinates are longitude-latitude decimal degrees, in which case distances are measured in kilometers; if coords is a SpatialPoints object, the value is taken from the object itself}
  \item{...}{additional arguments affecting the output produced}
  \item{x}{an object of class \code{nb}}
  \item{zero.policy}{default NULL, use global option value; if FALSE stop with error for any empty neighbour sets}
  \item{scale}{passed through to \code{stem()} for control of plot length}
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{plot.nb}}}

\examples{
if (require(rgdal, quietly=TRUE)) {
example(columbus, package="spData")
coords <- coordinates(columbus)
col.gal.nb
summary(col.gal.nb, coords)
col.listw <- nb2listw(col.gal.nb, style="W")
col.listw
summary(col.listw)
}}
\keyword{spatial}
