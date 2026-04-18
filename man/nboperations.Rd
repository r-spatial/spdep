% Copyright 2001, 2025 by Nicholas Lewin-Koh and Roger Bivand
\name{nb.set.operations}
\alias{intersect.nb}
\alias{union.nb}
\alias{setdiff.nb}
\alias{complement.nb}
\title{Set operations on neighborhood objects}
\description{
 Performs set operations on neighbors list objects union, intersection, (asymmetric!) difference and complement; note that from 1.3-9, \code{setdiff.nb} behaves like \code{\link[base]{setdiff}}, where the ordering of the objects being compared matters - (asymmetric!) difference.
}
\usage{
intersect.nb(nb.obj1,nb.obj2)
union.nb(nb.obj1,nb.obj2)
setdiff.nb(nb.obj1,nb.obj2)
complement.nb(nb.obj)
}

\arguments{
  \item{nb.obj}{a neighbor list created from any of the neighborhood
    list funtions}
  \item{nb.obj1}{a neighbor list created from any of the neighborhood
    list funtions}
  \item{nb.obj2}{a neighbor list created from any of the neighborhood
    list funtions}
}
\details{
  These functions perform set operations on each element of a
  neighborlist. The arguments must be neighbor lists created from the
  same spatial objects, checked by the region.id attributes being required
  to be identical. 
}
\value{
A new neighborlist created from the set operations on the input neighbor list(s)
}

\author{Nicholas Lewin-Koh \email{nikko@hailmail.net}}

\seealso{\code{\link{intersect.nb}}, \code{\link{union.nb}},
\code{\link{setdiff.nb}}}

\examples{
columbus <- st_read(system.file("shapes/columbus.gpkg", package="spData")[1], quiet=TRUE)
col.gal.nb <- read.gal(system.file("weights/columbus.gal", package="spData")[1])
coords <- st_coordinates(st_centroid(columbus))
col.tri.nb <- tri2nb(coords)
oldpar <- par(mfrow=c(1,2))
if (require("dbscan", quietly=TRUE)) {
  col.soi.nb <- graph2nb(soi.graph(col.tri.nb, coords))
  plot(st_geometry(columbus), border="grey")
  plot(col.soi.nb, coords, add=TRUE)
  title(main="Sphere of Influence Graph", cex.main=0.7)
  plot(st_geometry(columbus), border="grey")
  plot(complement.nb(col.soi.nb), coords, add=TRUE)
  title(main="Complement of Sphere of Influence Graph", cex.main=0.7)
}
par(mfrow=c(2,2))
col2 <- droplinks(col.gal.nb, 21)
plot(intersect.nb(col.gal.nb, col2), coords)
title(main="Intersect", cex.main=0.7)
plot(union.nb(col.gal.nb, col2), coords)
title(main="Union", cex.main=0.7)
plot(setdiff.nb(col.gal.nb, col2), coords)
title(main="Set diff", cex.main=0.7)
par(oldpar)
}
\keyword{spatial}
