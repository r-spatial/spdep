% Copyright 2001-2021 by Roger S. Bivand
\name{dnearneigh}
\alias{dnearneigh}
\title{Neighbourhood contiguity by distance}
\description{
The function identifies neighbours of region points by Euclidean distance between lower (greater than or equal to (changed from version 1.1-7)) and upper (less than or equal to) bounds, or with longlat = TRUE, by Great Circle distance in kilometers.
}
\usage{
dnearneigh(x, d1, d2, row.names = NULL, longlat = NULL, bounds=c("GE", "LE"),
 prefer_heuristic=TRUE, heuristic_n=1000, heuristic_fuzz=0.02, rtree=FALSE,
 exact_heuristic=FALSE, symtest=FALSE)
}
\arguments{
  \item{x}{matrix of point coordinates or a SpatialPoints object or an \code{sf} or \code{sfc} points object}
  \item{d1}{lower distance bound}
  \item{d2}{upper distance bound}
  \item{row.names}{character vector of region ids to be added to the neighbours list as attribute \code{region.id}, default \code{seq(1, nrow(x))}}
  \item{longlat}{TRUE if point coordinates are longitude-latitude decimal degrees, in which case distances are measured in kilometers; if x is a SpatialPoints object, the value is taken from the object itself, and overrides this argument if not NULL}
  \item{bounds}{character vector of length 2, default \code{c("GE", "LE")}, the first element may also be \code{"GE"}, the second \code{"LT"}; the first bound default was changed from \code{"GT"} to \code{"GE"} in release 1.1-7}
  \item{prefer_heuristic}{default TRUE, if TRUE, use \pkg{sf} \code{st_buffer()} and \code{st_intersects()} to find neighbours, or if \code{rtree=} is TRUE, uses \pkg{rtree} once submitted to CRAN.}
  \item{heuristic_n}{default 1000, do not use heuristics even if requested below this value.}
  \item{heuristic_fuzz}{default 0.02, for \code{st_buffer()} band taken as \code{d1} (if not 0) and \code{d2} times 1 +/- \code{heuristic_fuzz} because buffer imprecision increases with buffer radius for fixed \code{nQuadSegs} here always default 30, if \pkg{rtree} an absolute value in the metric of the points because that is how the R-tree is constructed; not used if \code{exact_heuristic} is FALSE.}
  \item{rtree}{default FALSE, may be used if \pkg{rtree} (\url{https://github.com/akoyabio/rtree}, \code{@BH_1.75} branch is installed). Code in CRAN package commented out, available on Github in branch \code{@rtree}.}
  \item{exact_heuristic}{default TRUE, if FALSE, use the raw buffer or rtree searcch without correcting for approximation.}
  \item{symtest}{Default FALSE; before release 1.1-7, TRUE - run symmetry check on output object, costly with large numbers of points.}
}
\value{
The function returns a list of integer vectors giving the region id numbers
for neighbours satisfying the distance criteria. See \code{\link{card}} for details of \dQuote{nb} objects.
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{knearneigh}}, \code{\link{card}}}

\examples{
columbus <- st_read(system.file("shapes/columbus.shp", package="spData")[1], quiet=TRUE)
coords <- st_centroid(st_geometry(columbus), of_largest_polygon=TRUE)
rn <- row.names(columbus)
k1 <- knn2nb(knearneigh(coords))
all.linked <- max(unlist(nbdists(k1, coords)))
col.nb.0.all <- dnearneigh(coords, 0, all.linked, row.names=rn)
summary(col.nb.0.all, coords)
opar <- par(no.readonly=TRUE)
plot(st_geometry(columbus), border="grey", reset=FALSE,
 main=paste("Distance based neighbours 0-",  format(all.linked), sep=""))
plot(col.nb.0.all, coords, add=TRUE)
par(opar)
(sfc_obj <- st_centroid(st_geometry(columbus)))
col.nb.0.all_sf <- dnearneigh(sfc_obj, 0, all.linked, row.names=rn)
all.equal(col.nb.0.all, col.nb.0.all_sf, check.attributes=FALSE)
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
title(main="Differences Euclidean/Great Circle")

xy1 <- SpatialPoints((as.data.frame(state.center))[m50.48,],
  proj4string=CRS("+proj=longlat +ellps=GRS80"))
gck1a <- knn2nb(knearneigh(xy1, k=1))
all.linked <- max(unlist(nbdists(gck1a, xy1)))
gc.nb <- dnearneigh(xy1, 0, all.linked)
summary(gc.nb, xy1, scale=0.5)
}
\keyword{spatial}
