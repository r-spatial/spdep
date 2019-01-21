% Copyright 2001 by Roger S. Bivand
\name{tri2nb}
\alias{tri2nb}
\title{Neighbours list from tri object}
\description{
The function uses the \code{deldir} package to convert a matrix of
two-dimensional coordinates into a neighbours list of class \code{nb}
with a list of integer vectors containing neighbour region number
ids. 
}
\usage{
tri2nb(coords, row.names = NULL)
}
\arguments{
  \item{coords}{matrix of point coordinates with two columns}
  \item{row.names}{character vector of region ids to be added to the neighbours list as attribute \code{region.id}, default \code{seq(1, nrow(x))}}
}

\details{
If coordinates are duplicated, this function cannot be used.
If the coordinates are from a grid, then they need to be ordered
such that the first three are not collinear, so that the first triangle
can be constructed. This can be achieved by randomising the order of
the coordinates (possibly several times), and then re-ordering the order
of the data to match the new order of the neighbour list - if this fix
is used, remember to re-order the row.names argument as well as the
coordinates! Please also note that triangulation of grid points will
give arbitrary diagonal neighbours, which may not be a sensible outcome,
and \code{dnearneigh()} may serve better where \code{tri2nb()} cannot
be used.
}

\value{
The function returns an object of class \code{nb} with a list of integer vectors containing neighbour region number ids.
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{knn2nb}}, \code{\link{dnearneigh}}, \code{\link{cell2nb}}}

\examples{
run <- FALSE
if (require("sf", quietly=TRUE)) run <- TRUE
if (run) {
columbus <- st_read(system.file("shapes/columbus.shp", package="spData")[1], quiet=TRUE)
}
if (run) {
coords <- coordinates(as(columbus, "Spatial"))
ind <- row.names(as(columbus, "Spatial"))
col.tri.nb <- tri2nb(coords, row.names=ind)
plot(st_geometry(columbus), border="grey")
plot(col.tri.nb, coords, add=TRUE)
title(main="Raw triangulation links", cex.main=0.6)
}
if (run) {
x <- seq(0,1,0.1)
y <- seq(0,2,0.2)
xy <- expand.grid(x, y)
try(xy.nb <- tri2nb(xy))
seed <- 1234
xid <- sample(1:nrow(xy))
xy.nb <- tri2nb(xy[xid,])
plot(xy.nb, xy[xid,])
}}
\keyword{spatial}
