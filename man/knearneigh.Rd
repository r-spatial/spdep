% Copyright 2001-14 by Roger S. Bivand
\name{knearneigh}
\alias{knearneigh}
\title{K nearest neighbours for spatial weights}
\description{
The function returns a matrix with the indices of points belonging to the set of the k nearest neighbours of each other. If longlat = TRUE, Great Circle distances are used. A warning will be given if identical points are found.
}
\usage{
knearneigh(x, k=1, longlat = NULL, use_kd_tree=TRUE)
}

\arguments{
  \item{x}{matrix of point coordinates, an object inheriting from SpatialPoints or an \code{"sf"} or \code{"sfc"} object; if the \code{"sf"} or \code{"sfc"} object geometries are in geographical coordinates (\code{sf::st_is_longlat(x) == TRUE} and \code{sf::sf_use_s2() == TRUE}), \pkg{s2} will be used to find the neighbours because it uses spatial indexing \url{https://github.com/r-spatial/s2/issues/125} as opposed to the legacy method which uses brute-force}
  \item{k}{number of nearest neighbours to be returned}
  \item{longlat}{TRUE if point coordinates are longitude-latitude decimal degrees, in which case distances are measured in kilometers; if x is a SpatialPoints object, the value is taken from the object itself; longlat will override \code{kd_tree}}
  \item{use_kd_tree}{logical value, if the \pkg{dbscan} package is available, use for finding k nearest neighbours when longlat is FALSE, and when there are no identical points; from \url{https://github.com/r-spatial/spdep/issues/38}, the input data may have more than two columns if \pkg{dbscan} is used}
}
\details{
The underlying legacy C code is based on the \code{knn} function in the \pkg{class} package. 
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
\code{\link{knn2nb}}, \code{\link[dbscan]{kNN}}}

\examples{
columbus <- st_read(system.file("shapes/columbus.shp", package="spData")[1], quiet=TRUE)
coords <- st_centroid(st_geometry(columbus), of_largest_polygon=TRUE)
col.knn <- knearneigh(coords, k=4)
plot(st_geometry(columbus), border="grey")
plot(knn2nb(col.knn), coords, add=TRUE)
title(main="K nearest neighbours, k = 4")
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
summary(llk4.nb, xy, longlat=TRUE, scale=0.5)
summary(gck4.nb, xy, longlat=TRUE, scale=0.5)

xy1 <- SpatialPoints((as.data.frame(state.center))[m50.48,],
  proj4string=CRS("+proj=longlat +ellps=GRS80"))
gck4a.nb <- knn2nb(knearneigh(xy1, k=4))
summary(gck4a.nb, xy1, scale=0.5)

xy1 <- st_as_sf((as.data.frame(state.center))[m50.48,], coords=1:2,
  crs=st_crs("+proj=longlat +ellps=GRS80"))
old_use_s2 <- sf_use_s2()
sf_use_s2(TRUE)
system.time(gck4a.nb <- knn2nb(knearneigh(xy1, k=4)))
summary(gck4a.nb, xy1, scale=0.5)
sf_use_s2(FALSE)
system.time(gck4a.nb <- knn2nb(knearneigh(xy1, k=4)))
summary(gck4a.nb, xy1, scale=0.5)
sf_use_s2(old_use_s2)

# https://github.com/r-spatial/spdep/issues/38
if (require("dbscan", quietly=TRUE)) {
  set.seed(1)
  x <- cbind(runif(50), runif(50), runif(50))
  out <- knearneigh(x, k=5)
  knn2nb(out)
  try(out <- knearneigh(rbind(x, x[1:10,]), k=5))
}
}
\keyword{spatial}
