% Copyright 2001-3 by Roger S. Bivand
\name{diffnb}
\alias{diffnb}
\title{Differences between neighbours lists}
\description{
The function finds differences between lists of neighbours, returning a \code{nb} neighbour list of those found
}
\usage{
diffnb(x, y, verbose=NULL)
}
\arguments{
  \item{x}{an object of class \code{nb}}
  \item{y}{an object of class \code{nb}}
  \item{verbose}{default NULL, use global option value; report regions ids taken from object attribute "region.id" with differences}
}

\value{
A neighbours list with class \code{nb}
}

\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\examples{
run <- FALSE
if (require("sf", quietly=TRUE)) run <- TRUE
if (run) {
columbus <- st_read(system.file("shapes/columbus.shp", package="spData")[1], quiet=TRUE)
}
if (run) {
coords <- st_coordinates(st_centroid(st_geometry(columbus)))
rn <- row.names(columbus)
knn1 <- knearneigh(coords, 1)
knn2 <- knearneigh(coords, 2)
nb1 <- knn2nb(knn1, row.names=rn)
nb2 <- knn2nb(knn2, row.names=rn)
diffs <- diffnb(nb2, nb1)
opar <- par(no.readonly=TRUE)
plot(st_geometry(columbus), border="grey", reset=FALSE,
 main="Plot of first (black) and second (red)\nnearest neighbours")
plot(nb1, coords, add=TRUE)
plot(diffs, coords, add=TRUE, col="red", lty=2)
par(opar)
}
}
\keyword{spatial}
