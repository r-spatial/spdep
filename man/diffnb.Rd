% Copyright 2001-3, 25 by Roger S. Bivand
\name{diffnb}
\alias{diffnb}
\title{Differences between neighbours lists}
\description{
The function finds symmetric differences between lists of neighbours (ordering of objects does not matter), returning a \code{nb} neighbour list of those found
}
\usage{
diffnb(x, y, verbose=NULL, legacy=TRUE)
}
\arguments{
  \item{x}{an object of class \code{nb}}
  \item{y}{an object of class \code{nb}}
  \item{verbose}{default NULL, use global option value; report regions ids taken from object attribute "region.id" with differences, ignored when \code{legacy=} is FALSE}
  \item{legacy}{default TRUE; use legacy code; if FALSE use differences between sparse matrix representations}
}

\value{
A neighbours list with class \code{nb}
}

\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\examples{
columbus <- st_read(system.file("shapes/columbus.gpkg", package="spData")[1], quiet=TRUE)
coords <- st_centroid(st_geometry(columbus), of_largest_polygon=TRUE)
rn <- row.names(columbus)
knn1 <- knearneigh(coords, 1)
knn2 <- knearneigh(coords, 2)
nb1 <- knn2nb(knn1, row.names=rn)
nb2 <- knn2nb(knn2, row.names=rn)
diffs <- diffnb(nb2, nb1, legacy=TRUE)
diffsl <- diffnb(nb2, nb1, legacy=FALSE)
all.equal(diffs, diffsl)
# call attribute fifference expected
diffsd <- union.nb(setdiff.nb(nb1, nb2), setdiff.nb(nb2, nb1))
all.equal(diffs, diffsd)
# call attribute fifference expected
opar <- par(no.readonly=TRUE)
plot(st_geometry(columbus), border="grey", reset=FALSE,
 main="Plot of first (black) and second (red)\nnearest neighbours")
plot(nb1, coords, add=TRUE)
plot(diffs, coords, add=TRUE, col="red", lty=2)
par(opar)
}
\keyword{spatial}
