% Copyright 2001-2021 by Roger S. Bivand
\name{dnearneigh}
\alias{dnearneigh}
\title{Neighbourhood contiguity by distance}
\description{
The function identifies neighbours of region points by Euclidean distance between lower (greater than or equal to (changed from version 1.1-7)) and upper (less than or equal to) bounds, or with longlat = TRUE, by Great Circle distance in kilometers.
}
\usage{
dnearneigh(x, d1, d2, row.names = NULL, longlat = NULL, bounds=c("GE", "LE"),
 use_kd_tree=TRUE, symtest=FALSE, use_s2=FALSE, max_cells=200, dwithin=FALSE)
}
\arguments{
  \item{x}{matrix of point coordinates, an object inheriting from SpatialPoints or an \code{"sf"} or \code{"sfc"} object; if the \code{"sf"} or \code{"sfc"} object geometries are in geographical coordinates (\code{use_s2=FALSE}, \code{sf::st_is_longlat(x) == TRUE} and \code{sf::sf_use_s2() == TRUE}), \pkg{s2} will be used to find the neighbours because it will (we hope) use spatial indexing \url{https://github.com/r-spatial/s2/issues/125} as opposed to the legacy method which uses brute-force (at present \pkg{s2} also uses brute-force)}
  \item{d1}{lower distance bound}
  \item{d2}{upper distance bound}
  \item{row.names}{character vector of region ids to be added to the neighbours list as attribute \code{region.id}, default \code{seq(1, nrow(x))}}
  \item{longlat}{TRUE if point coordinates are longitude-latitude decimal degrees, in which case distances are measured in kilometers; if x is a SpatialPoints object, the value is taken from the object itself, and overrides this argument if not NULL}
  \item{bounds}{character vector of length 2, default \code{c("GE", "LE")}, (GE: greater than or equal to, LE: less than or equal to) that is the finite and closed interval \code{[d1, d2]}, \code{d1 <= x <= d2}. The first element may also be \code{"GT"} (GT: greater than), the second \code{"LT"} (LT: less than) for finite, open intervals excluding the bounds; the first bound default was changed from \code{"GT"} to \code{"GE"} in release 1.1-7. When creating multiple distance bands, finite, half-open right-closed intervals may be used until the final interval to avoid overlapping on bounds: \code{"GE", "LT"}, that is \code{[d1, d2)}, \code{d1 <= x < d2}}
  \item{use_kd_tree}{default TRUE, if TRUE, use \pkg{dbscan} \code{\link[dbscan]{frNN}} if available (permitting 3D distances).}
  \item{symtest}{Default FALSE; before release 1.1-7, TRUE - run symmetry check on output object, costly with large numbers of points.}
  \item{use_s2}{default=FALSE, as of \pkg{s2} 1.0-5, distance bound compuations do not use spatial indexing so even if \code{sf::sf_use_s2()} is TRUE, \pkg{s2} will not be used for distances on the sphere for \code{"sf"} or \code{"sfc"} objects; the default will change when/if spatial indexing becomes available. Until then, setting \code{set.coresOption} to a positive integer will experimentaly use that number of cores on \code{.Platform$OS.type== "unix"} systems.}
  \item{max_cells}{default 200; if \code{use_s2=TRUE}, \code{sf::st_is_longlat(x) == TRUE}, \code{sf::sf_use_s2() == TRUE} and \code{dwithin=FALSE}, set \code{max_cells} in \code{s2::s2_buffer_cells()}}
  \item{dwithin}{default FALSE to use \code{s2::s2_buffer_cells()} and \code{s2::s2_intersects_matrix()}, if TRUE, use \code{s2::s2_dwithin_matrix()}, both if \code{use_s2=TRUE}, \code{sf::st_is_longlat(x) == TRUE} and \code{sf::sf_use_s2() == TRUE}; both yield the same lists of neighbours.}
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
(all.linked <- max(unlist(nbdists(llk1, xy, longlat=FALSE))))
ll.nb <- dnearneigh(xy, 0, all.linked, longlat=FALSE)
summary(ll.nb, xy, longlat=TRUE, scale=0.5)
gck1 <- knn2nb(knearneigh(xy, k=1, longlat=TRUE))
(all.linked <- max(unlist(nbdists(gck1, xy, longlat=TRUE))))
gc.nb <- dnearneigh(xy, 0, all.linked, longlat=TRUE)
summary(gc.nb, xy, longlat=TRUE, scale=0.5)
plot(ll.nb, xy)
plot(diffnb(ll.nb, gc.nb), xy, add=TRUE, col="red", lty=2)
title(main="Differences Euclidean/Great Circle")

xy1 <- SpatialPoints((as.data.frame(state.center))[m50.48,],
  proj4string=CRS("+proj=longlat +ellps=GRS80"))
gck1a <- knn2nb(knearneigh(xy1, k=1))
(all.linked <- max(unlist(nbdists(gck1a, xy1))))
gc.nb <- dnearneigh(xy1, 0, all.linked)
summary(gc.nb, xy1, scale=0.5)

xy1 <- st_as_sf((as.data.frame(state.center))[m50.48,], coords=1:2,
  crs=st_crs("+proj=longlat +ellps=GRS80"))
old_use_s2 <- sf_use_s2()
sf_use_s2(TRUE)
gck1b <- knn2nb(knearneigh(xy1, k=1))
system.time(o <- nbdists(gck1b, xy1))
(all.linked <- max(unlist(o)))
# use s2 brute-force buffer/intersect approach (with two passes if d1 > 0)
system.time(gc.nb <- dnearneigh(xy1, 0, all.linked, use_s2=TRUE))
summary(gc.nb, xy1, scale=0.5)
# use s2 brute-force buffer/intersect approach with smaller max_cells
system.time(gc.nb <- dnearneigh(xy1, 0, all.linked, use_s2=TRUE, max_cells=500))
summary(gc.nb, xy1, scale=0.5)
# use s2 brute-force dwithin_matrix approach 
system.time(gc.nb <- dnearneigh(xy1, 0, all.linked, use_s2=TRUE, dwithin=TRUE))
summary(gc.nb, xy1, scale=0.5)
# use legacy symmetric brute-force approach
system.time(gc.nb <- dnearneigh(xy1, 0, all.linked))
summary(gc.nb, xy1, scale=0.5)
sf_use_s2(old_use_s2)
}
\keyword{spatial}
