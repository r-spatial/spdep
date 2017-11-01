% Copyright 2001-10 by Roger S. Bivand
\name{poly2nb}
\alias{poly2nb}
\title{Construct neighbours list from polygon list}
\description{
The function builds a neighbours list based on regions with contiguous boundaries, that is sharing one or more boundary point. The current function is in part interpreted and may run slowly for many regions or detailed boundaries, but from 0.2-16 should not fail because of lack of memory when single polygons are built of very many border coordinates.}
\usage{
poly2nb(pl, row.names = NULL, snap=sqrt(.Machine$double.eps),
 queen=TRUE, useC=TRUE, foundInBox=NULL)
}

\arguments{
  \item{pl}{list of polygons of class extending \code{SpatialPolygons}}
  \item{row.names}{character vector of region ids to be added to the neighbours list as attribute \code{region.id}, default \code{seq(1, nrow(x))}; if \code{polys} has a \code{region.id} attribute, it is copied to the neighbours list.}
  \item{snap}{boundary points less than \code{snap} distance apart are considered to indicate contiguity}
  \item{queen}{if TRUE, a single shared boundary point meets the contiguity condition, if FALSE, more than one shared point is required; note that more than one shared boundary point does not necessarily mean a shared boundary line}
  \item{useC}{default TRUE, doing the work loop in C, may be set to false to revert to R code calling two C functions in an \code{n*k} work loop, where \code{k} is the average number of candidate neighbours}
  \item{foundInBox}{default NULL using R code, possibly parallelised if a \pkg{snow} cluster is available, otherwise a list of length \code{(n-1)} with integer vectors of candidate neighbours \code{(j > i)}, or NULL if all candidates were \code{(j < i)} (as created by the \code{poly_findInBoxGEOS} function in \pkg{rgeos} for clean polygons)}
}
\value{
A neighbours list with class \code{nb}. See \code{\link{card}} for details of \dQuote{nb} objects.
}
\note{
From 0.5-8, the function includes faster bounding box indexing and other improvements contributed by Micah Altman. If a cluster is provided using \code{set.ClusterOption}, it will be used for finding candidate bounding box overlaps for exact testing for contiguity.
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no} with contributions from Micah Altman}

\seealso{\code{\link{summary.nb}}, \code{\link{card}}}

\examples{
if (require(rgdal, quietly=TRUE)) {
example(columbus, package="spData")
coords <- coordinates(columbus)
xx <- poly2nb(columbus)
dxx <- diffnb(xx, col.gal.nb)
plot(columbus, border="grey")
plot(col.gal.nb, coords, add=TRUE)
plot(dxx, coords, add=TRUE, col="red")
title(main=paste("Differences (red) in Columbus GAL weights (black)",
 "and polygon generated queen weights", sep="\n"))
xxx <- poly2nb(columbus, queen=FALSE)
dxxx <- diffnb(xxx, col.gal.nb)
plot(columbus, border = "grey")
plot(col.gal.nb, coords, add = TRUE)
plot(dxxx, coords, add = TRUE, col = "red")
title(main=paste("Differences (red) in Columbus GAL weights (black)",
 "and polygon generated rook weights", sep="\n"))
cards <- card(xx)
maxconts <- which(cards == max(cards))
if(length(maxconts) > 1) maxconts <- maxconts[1]
fg <- rep("grey", length(cards))
fg[maxconts] <- "red"
fg[xx[[maxconts]]] <- "green"
plot(columbus, col=fg)
title(main="Region with largest number of contiguities")
example(nc.sids, package="spData")
system.time(xxnb <- poly2nb(nc.sids))
plot(nc.sids)
plot(xxnb, coordinates(nc.sids), add=TRUE, col="blue")
}
}
\keyword{spatial}
