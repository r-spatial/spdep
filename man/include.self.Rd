% Copyright 2001 by Roger S. Bivand
\name{include.self}
\alias{include.self}
\alias{remove.self}
\title{Include self in neighbours list}
\description{
The function includes the region itself in its own list of neighbours, and sets attribute "self.included" to TRUE; \code{remove.self} reverts the effects of \code{include.self}.
}
\usage{
include.self(nb)
remove.self(nb)
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
columbus <- st_read(system.file("shapes/columbus.shp", package="spData")[1], quiet=TRUE)
col.gal.nb <- read.gal(system.file("weights/columbus.gal", package="spData")[1])
coords <- st_coordinates(st_centroid(columbus))
summary(col.gal.nb, coords)
summary(include.self(col.gal.nb), coords)
summary(remove.self(include.self(col.gal.nb)), coords)
}
\keyword{spatial}
