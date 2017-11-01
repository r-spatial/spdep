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
if (require(rgdal, quietly=TRUE)) {
example(columbus, package="spData")
coords <- coordinates(columbus)
summary(col.gal.nb, coords)
summary(include.self(col.gal.nb), coords)
}
}
\keyword{spatial}
