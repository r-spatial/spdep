% Copyright 2001-14 by Roger S. Bivand
\name{knearneigh}
\alias{knearneigh}
\title{K nearest neighbours for spatial weights}
\description{
The function returns a matrix with the indices of points belonging to the set of the k nearest neighbours of each other. If longlat = TRUE, Great Circle distances are used. A warning will be given if identical points are found.
}
\usage{
knearneigh(x, k=1, longlat = NULL, RANN=TRUE)
}

\arguments{
  \item{x}{matrix of point coordinates or a SpatialPoints object}
  \item{k}{number of nearest neighbours to be returned}
  \item{longlat}{TRUE if point coordinates are longitude-latitude decimal degrees, in which case distances are measured in kilometers; if x is a SpatialPoints object, the value is taken from the object itself}
  \item{RANN}{logical value, if the RANN package is available, use for finding k nearest neighbours when longlat is FALSE, and when there are no identical points}
}
\details{
The underlying C code is based on the \code{knn} function in the \pkg{class} package.
}
\value{
A list of class \code{knn}
  \item{nn}{integer matrix of region number ids}
  \item{np}{number of input points}
  \item{k}{input required k}
  \item{dimension}{number of columns of x}
  \item{x}{input coordinates}
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link[class]{knn}}, \code{\link{dnearneigh}},
\code{\link{knn2nb}}, \code{\link[RANN]{nn2}}}

\examples{
if (require(rgdal, quietly=TRUE)) {
example(columbus, package="spData")
coords <- coordinates(columbus)
col.knn <- knearneigh(coords, k=4)
plot(columbus, border="grey")
plot(knn2nb(col.knn), coords, add=TRUE)
title(main="K nearest neighbours, k = 4")
}
data(state)
us48.fipsno <- read.geoda(system.file("etc/weights/us48.txt",
 package="spdep")[1])
if (as.numeric(paste(version$major, version$minor, sep="")) < 19) {
 m50.48 <- match(us48.fipsno$"State.name", state.name)
} else {
 m50.48 <- match(us48.fipsno$"State_name", state.name)
}
xy <- as.matrix(as.data.frame(state.center))[m50.48,]
llk4.nb <- knn2nb(knearneigh(xy, k=4, longlat=FALSE))
gck4.nb <- knn2nb(knearneigh(xy, k=4, longlat=TRUE))
plot(llk4.nb, xy)
plot(diffnb(llk4.nb, gck4.nb), xy, add=TRUE, col="red", lty=2)
title(main="Differences between Euclidean and Great Circle k=4 neighbours")
summary(llk4.nb, xy, longlat=TRUE)
summary(gck4.nb, xy, longlat=TRUE)

xy1 <- SpatialPoints((as.data.frame(state.center))[m50.48,],
  proj4string=CRS("+proj=longlat +ellps=GRS80"))
gck4a.nb <- knn2nb(knearneigh(xy1, k=4))
summary(gck4a.nb, xy1)
}
\keyword{spatial}
