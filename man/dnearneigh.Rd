% Copyright 2001-2014 by Roger S. Bivand
\name{dnearneigh}
\alias{dnearneigh}
\title{Neighbourhood contiguity by distance}
\description{
The function identifies neighbours of region points by Euclidean distance between lower (greater than) and upper (less than or equal to) bounds, or with longlat = TRUE, by Great Circle distance in kilometers.
}
\usage{
dnearneigh(x, d1, d2, row.names = NULL, longlat = NULL, bounds=c("GT", "LE"))
}
\arguments{
  \item{x}{matrix of point coordinates or a SpatialPoints object}
  \item{d1}{lower distance bound}
  \item{d2}{upper distance bound}
  \item{row.names}{character vector of region ids to be added to the neighbours list as attribute \code{region.id}, default \code{seq(1, nrow(x))}}
  \item{longlat}{TRUE if point coordinates are longitude-latitude decimal degrees, in which case distances are measured in kilometers; if x is a SpatialPoints object, the value is taken from the object itself, and overrides this argument if not NULL}
  \item{bounds}{character vector of length 2, default \code{c("GT", "LE")}, the first element may also be \code{"GE"}, the second \code{"LT"}}
}
\value{
The function returns a list of integer vectors giving the region id numbers
for neighbours satisfying the distance criteria. See \code{\link{card}} for details of \dQuote{nb} objects.
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{knearneigh}}, \code{\link{card}}}

\examples{
if (require(rgdal, quietly=TRUE)) {
example(columbus, package="spData")
coords <- coordinates(columbus)
rn <- sapply(slot(columbus, "polygons"), function(x) slot(x, "ID"))
k1 <- knn2nb(knearneigh(coords))
all.linked <- max(unlist(nbdists(k1, coords)))
col.nb.0.all <- dnearneigh(coords, 0, all.linked, row.names=rn)
summary(col.nb.0.all, coords)
plot(columbus, border="grey")
plot(col.nb.0.all, coords, add=TRUE)
title(main=paste("Distance based neighbours 0-",  format(all.linked),
 " distance units", sep=""))
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
llk1 <- knn2nb(knearneigh(xy, k=1, longlat=FALSE))
all.linked <- max(unlist(nbdists(llk1, xy, longlat=FALSE)))
ll.nb <- dnearneigh(xy, 0, all.linked, longlat=FALSE)
summary(ll.nb, xy, longlat=TRUE, scale=0.5)
gck1 <- knn2nb(knearneigh(xy, k=1, longlat=TRUE))
all.linked <- max(unlist(nbdists(gck1, xy, longlat=TRUE)))
gc.nb <- dnearneigh(xy, 0, all.linked, longlat=TRUE)
summary(gc.nb, xy, longlat=TRUE, scale=0.5)
plot(ll.nb, xy)
plot(diffnb(ll.nb, gc.nb), xy, add=TRUE, col="red", lty=2)
title(main="Differences between Euclidean and Great Circle neighbours")

xy1 <- SpatialPoints((as.data.frame(state.center))[m50.48,],
  proj4string=CRS("+proj=longlat +ellps=GRS80"))
gck1a <- knn2nb(knearneigh(xy1, k=1))
all.linked <- max(unlist(nbdists(gck1a, xy1)))
gc.nb <- dnearneigh(xy1, 0, all.linked)
summary(gc.nb, xy1, scale=0.5)
}
\keyword{spatial}
