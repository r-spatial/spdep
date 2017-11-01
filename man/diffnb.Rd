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
if (require(rgdal, quietly=TRUE)) {
example(columbus, package="spData")
coords <- coordinates(columbus)
rn <- sapply(slot(columbus, "polygons"), function(x) slot(x, "ID"))
knn1 <- knearneigh(coords, 1)
knn2 <- knearneigh(coords, 2)
nb1 <- knn2nb(knn1, row.names=rn)
nb2 <- knn2nb(knn2, row.names=rn)
diffs <- diffnb(nb2, nb1)
plot(columbus, border="grey")
plot(nb1, coords, add=TRUE)
plot(diffs, coords, add=TRUE, col="red", lty=2)
title(main="Plot of first (black) and second (red)\nnearest neighbours")
}
}
\keyword{spatial}
