% Copyright 2001 by Roger S. Bivand
\name{include.self}
\alias{include.self}
\title{Include self in neighbours list}
\description{
The function includes the region itself in its own list of neighbours, and sets attribute "self.included" to TRUE.
}
\usage{
include.self(nb)
}
\arguments{
  \item{nb}{input neighbours list of class \code{nb}}
}
\value{
The function returns an object of class \code{nb} with a list of integer vectors containing neighbour region number ids; attribute "self.included" is set to TRUE.
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{summary.nb}}}

\examples{
run <- FALSE
if (require("sf", quietly=TRUE)) run <- TRUE
if (run) {
columbus <- st_read(system.file("shapes/columbus.shp", package="spData")[1], quiet=TRUE)
col.gal.nb <- read.gal(system.file("weights/columbus.gal", package="spData")[1])
}
if (run) {
coords <- coordinates(as(columbus, "Spatial"))
summary(col.gal.nb, coords)
}
if (run) {
summary(include.self(col.gal.nb), coords)
}
}
\keyword{spatial}
