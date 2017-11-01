% Copyright 2001 by Nicholas Lewin-Koh
\name{nb.set.operations}
\alias{intersect.nb}
\alias{union.nb}
\alias{setdiff.nb}
\alias{complement.nb}
\title{Set operations on neighborhood objects}
\description{
 Set operations on neighbors list objects
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
  same coordinates, and the region.id attributes must be identical. 
}
\value{
 \item{nb.obj}{A new neighborlist created from the set operations on the
               input neighbor list(s)} 
}

\author{Nicholas Lewin-Koh \email{nikko@hailmail.net}}

\seealso{\code{\link{intersect.nb}}, \code{\link{union.nb}},
\code{\link{setdiff.nb}}}

\examples{
if (require(rgdal, quietly=TRUE)) {
example(columbus, package="spData")
coords <- coordinates(columbus)
col.tri.nb <- tri2nb(coords)
oldpar <- par(mfrow=c(1,2))
col.soi.nb <- graph2nb(soi.graph(col.tri.nb, coords))
plot(columbus, border="grey")
plot(col.soi.nb, coords, add=TRUE)
title(main="Sphere of Influence Graph")
plot(columbus, border="grey")
plot(complement.nb(col.soi.nb), coords, add=TRUE)
title(main="Complement of Sphere of Influence Graph")
par(mfrow=c(2,2))
col2 <- droplinks(col.gal.nb, 21)
plot(intersect.nb(col.gal.nb, col2), coords)
title(main="Intersect")
plot(union.nb(col.gal.nb, col2), coords)
title(main="Union")
plot(setdiff.nb(col.gal.nb, col2), coords)
title(main="Set diff")
par(oldpar)
}
}
\keyword{spatial}
