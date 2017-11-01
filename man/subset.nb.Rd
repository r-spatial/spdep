% Copyright 2001 by Roger S. Bivand
\name{subset.nb}
\alias{subset.nb}
\title{Subset a neighbours list}
\description{
  The function subsets a neighbors list, retaining objects for which the subset argument vector is TRUE.
}
\usage{
\method{subset}{nb}(x, subset, ...)
}

\arguments{
  \item{x}{an object of class \code{nb}}
  \item{subset}{logical expression}
  \item{...}{generic function pass-through}
}
\value{
  The function returns an object of class \code{nb} with a list of
   integer vectors containing neighbour region number ids (compacted to
   run from 1:number of regions in subset).
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}
\seealso{\code{\link{nb2listw}}}

\examples{
if (require(rgdal, quietly=TRUE)) {
example(columbus, package="spData")
coords <- coordinates(columbus)
plot(col.gal.nb, coords)
to.be.dropped <- c(31, 34, 36, 39, 42, 46)
text(coords[to.be.dropped,1], coords[to.be.dropped,2], labels=to.be.dropped,
  pos=2, offset=0.3)
sub.col.gal.nb <- subset(col.gal.nb,
  !(1:length(col.gal.nb) \%in\% to.be.dropped))
plot(sub.col.gal.nb, coords[-to.be.dropped,], col="red", add=TRUE)
which(!(attr(col.gal.nb, "region.id") \%in\%
  attr(sub.col.gal.nb, "region.id")))
}}
\keyword{spatial}
