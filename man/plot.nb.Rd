% Copyright 2001-5 by Roger S. Bivand
\name{plot.nb}
\alias{plot.nb}
\alias{plot.listw}
\title{Plot a neighbours list}
\description{
A function to plot a neighbours list given point coordinates to represent the region in two dimensions; \code{plot.listw} is a wrapper that passes its neighbours component to \code{plot.nb}.
}
\usage{
\method{plot}{nb}(x, coords, col="black", points=TRUE, add=FALSE, arrows=FALSE,
 length=0.1, xlim=NULL, ylim=NULL, ...)
\method{plot}{listw}(x, coords, col="black", points=TRUE, add=FALSE, arrows=FALSE,
 length=0.1, xlim=NULL, ylim=NULL, ...)
}

\arguments{
  \item{x}{an object of class \code{nb} or (for \code{plot.listw}) class \code{listw}}
  \item{coords}{matrix of region point coordinates}
  \item{col}{plotting colour}
  \item{points}{(logical) add points to plot}
  \item{add}{(logical) add to existing plot}
  \item{arrows}{(logical) draw arrowheads for asymmetric neighbours}
  \item{length}{length in plot inches of arrow heads drawn for asymmetric
neighbours lists}
  \item{xlim, ylim}{plot window bounds}
  \item{...}{further graphical parameters as in \code{par(..)}}
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{summary.nb}}}

\examples{
if (require(rgdal, quietly=TRUE)) {
example(columbus, package="spData")
coords <- coordinates(columbus)
plot(col.gal.nb, coords)
title(main="GAL order 1 links with first nearest neighbours in red")
col.knn <- knearneigh(coords, k=1)
plot(knn2nb(col.knn), coords, add=TRUE, col="red", length=0.08)
}}
\keyword{spatial}
