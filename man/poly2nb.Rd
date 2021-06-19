% Copyright 2001-19 by Roger S. Bivand
\name{poly2nb}
\alias{poly2nb}
\title{Construct neighbours list from polygon list}
\description{
The function builds a neighbours list based on regions with contiguous boundaries, that is sharing one or more boundary point. The current function is in part interpreted and may run slowly for many regions or detailed boundaries, but from 0.2-16 should not fail because of lack of memory when single polygons are built of very many border coordinates.}
\usage{
poly2nb(pl, row.names = NULL, snap=sqrt(.Machine$double.eps),
 queen=TRUE, useC=TRUE, foundInBox=NULL, small_n=500, use_s2=FALSE,
 dfMaxLength=0)
}

\arguments{
  \item{pl}{list of polygons of class extending \code{SpatialPolygons}, or an \code{sf} or \code{sfc} object containing non-empty (multi-)polygon objects}
  \item{row.names}{character vector of region ids to be added to the neighbours list as attribute \code{region.id}, default \code{seq(1, nrow(x))}; if \code{pl} has \code{row.names}, they are used instead of the default sequence.}
  \item{snap}{boundary points less than \code{snap} distance apart are considered to indicate contiguity}
  \item{queen}{if TRUE, a single shared boundary point meets the contiguity condition, if FALSE, more than one shared point is required; note that more than one shared boundary point does not necessarily mean a shared boundary line}
  \item{useC}{default TRUE, doing the work loop in C, may be set to false to revert to R code calling two C functions in an \code{n*k} work loop, where \code{k} is the average number of candidate neighbours}
  \item{foundInBox}{default NULL using R code, possibly parallelised if a \pkg{snow} cluster is available, otherwise a list of length \code{(n-1)} with integer vectors of candidate neighbours \code{(j > i)}, or NULL if all candidates were \code{(j < i)} (as created by the \code{poly_findInBoxGEOS} function in \pkg{rgeos} for clean polygons)}
  \item{small_n}{default 500; if \code{foundInBox=} is not given, generate internally using \pkg{sf} code if the number of geometries \code{>= small_n}. If the geometries are not clean, set to a value larger than the geometry count.}
  \item{use_s2}{default FALSE, do not re-assign the coordinate reference system to the geometry bounding boxes, thus using planar bounding box intersection; if TRUE, use (currently fallible) spherical bounding box intersection.}
  \item{dfMaxLength}{default 0 metres, if a positive length in metres and \code{use_s2=TRUE}, use in \code{st_segmentize()} to add points to bounding box edges.}
}
\value{
A neighbours list with class \code{nb}. See \code{\link{card}} for details of \dQuote{nb} objects.
}
\note{
From 0.5-8, the function includes faster bounding box indexing and other improvements contributed by Micah Altman. If a cluster is provided using \code{set.ClusterOption}, it will be used for finding candidate bounding box overlaps for exact testing for contiguity.

Until 1.1-7, \code{sf} polygons included both start and end points, so could erroneously report rook neighbourhood where only queen was present, see \url{https://github.com/r-spatial/spdep/issues/50}.

From 1.1-9 with \pkg{sf} 1.0-0, \pkg{s2} is used in bounding box indexing internally when \code{pl} is in geographical coordinates. Because the topology engine of \pkg{s2} differs from the use of GEOS for planar coordinates by \pkg{sf}, some output differences may be expected. When \pkg{s2} is used, (\code{use_s2=TRUE} and \code{dfMaxLength > 0}) segmenting is done to add points to bounding box edges, making them more likely to intersect on the sphere, but currently increasing timing and memory use substantially. If segmenting not is done, but \code{use_s2=TRUE}, candidate neighbours may be missed and others exchanged. Since treating spherical geometries as planar is also questionable, it is not clear whether spherical contiguous polygon neighbours should simply follow neighbours found by treating the geometries as planar \url{https://github.com/r-spatial/s2/issues/125#issuecomment-864403372}. However, current advice is not to use \pkg{s2} for finding contiguity neighbours.
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no} with contributions from Micah Altman}

\seealso{\code{\link{summary.nb}}, \code{\link{card}}}

\examples{
columbus <- st_read(system.file("shapes/columbus.shp", package="spData")[1], quiet=TRUE)
col.gal.nb <- read.gal(system.file("weights/columbus.gal", package="spData")[1])
coords <- st_coordinates(st_centroid(st_geometry(columbus)))
xx <- poly2nb(as(columbus, "Spatial"))
dxx <- diffnb(xx, col.gal.nb)
plot(st_geometry(columbus), border="grey")
plot(col.gal.nb, coords, add=TRUE)
plot(dxx, coords, add=TRUE, col="red")
title(main=paste("Differences (red) in Columbus GAL weights (black)",
 "and polygon generated queen weights", sep="\n"), cex.main=0.6)
# poly2nb with sf sfc_MULTIPOLYGON objects
sf_xx <- poly2nb(columbus)
diffnb(sf_xx, xx)
sfc_xx <- poly2nb(st_geometry(columbus))
diffnb(sfc_xx, xx)
xxx <- poly2nb(as(columbus, "Spatial"), queen=FALSE)
dxxx <- diffnb(xxx, col.gal.nb)
plot(st_geometry(columbus), border = "grey")
plot(col.gal.nb, coords, add = TRUE)
plot(dxxx, coords, add = TRUE, col = "red")
title(main=paste("Differences (red) in Columbus GAL weights (black)",
 "and polygon generated rook weights", sep="\n"), cex.main=0.6)
cards <- card(xx)
maxconts <- which(cards == max(cards))
if(length(maxconts) > 1) maxconts <- maxconts[1]
fg <- rep("grey", length(cards))
fg[maxconts] <- "red"
fg[xx[[maxconts]]] <- "green"
plot(st_geometry(columbus), col=fg)
title(main="Region with largest number of contiguities", cex.main=0.6)
nc.sids <- st_read(system.file("shapes/sids.shp", package="spData")[1], quiet=TRUE)
system.time(xxnb <- poly2nb(nc.sids))
system.time(xxnb <- poly2nb(as(nc.sids, "Spatial")))
plot(st_geometry(nc.sids))
plot(xxnb, coordinates(as(nc.sids, "Spatial")), add=TRUE, col="blue")
sq <- st_polygon(list(rbind(c(0,0), c(1,0), c(1,1), c(0,1), c(0,0))))
sq2 <- sq + c(0,1)
sq3 <- sq + c(1,0)
sq4 <- sq + c(1,1)
gm <- st_sfc(list(sq, sq2, sq3, sq4))
df <- st_as_sf(gm, id=1:4)
plot(st_geometry(df))
text(st_coordinates(st_centroid(gm)), as.character(df$id))
unclass(poly2nb(df, queen = FALSE))
}
\keyword{spatial}
