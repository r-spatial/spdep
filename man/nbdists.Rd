% Copyright 2001 by Roger S. Bivand
\name{nbdists}
\alias{nbdists}
\title{Spatial link distance measures}
\description{
  Given a list of spatial neighbour links (a neighbours list of object type
  \code{nb}), the function returns the Euclidean distances along the links
  in a list of the same form as the neighbours list. If longlat = TRUE, Great
  Circle distances are used.
}
\usage{
nbdists(nb, coords, longlat = NULL)
}
\arguments{
  \item{nb}{an object of class \code{nb}}
  \item{coords}{matrix of point coordinates or a SpatialPoints object}
  \item{longlat}{TRUE if point coordinates are longitude-latitude decimal degrees, in which case distances are measured in kilometers; if coords is a SpatialPoints object, the value is taken from the object itself}
}
\value{
A list with class \code{nbdist}
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{summary.nb}}, \code{\link{nb2listw}}}

\examples{
run <- FALSE
if (require("sf", quietly=TRUE)) run <- TRUE
if (run) {
columbus <- st_read(system.file("shapes/columbus.shp", package="spData")[1], quiet=TRUE)
col.gal.nb <- read.gal(system.file("weights/columbus.gal", package="spData")[1])
coords <- coordinates(as(columbus, "Spatial"))
dlist <- nbdists(col.gal.nb, coords)
dlist <- lapply(dlist, function(x) 1/x)
stem(unlist(dlist))
}
}
\keyword{spatial}
