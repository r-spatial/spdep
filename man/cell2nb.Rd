% Copyright 2001-2021 by Roger S. Bivand
\name{cell2nb}
\alias{cell2nb}
\alias{vi2mrc}
\title{Generate neighbours list for grid cells}
\description{
 The function generates a list of neighbours for a grid of cells. Helper
functions are used to convert to and from the vector indices for row and 
column grid positions, and rook (shared edge) or queen (shared edge or 
vertex) neighbour definitions are applied by type. If torus is TRUE, the 
grid is mapped onto a torus, removing edge effects.
}
\usage{
cell2nb(nrow, ncol, type="rook", torus=FALSE, legacy=FALSE, x=NULL)
vi2mrc(i, nrow, ncol)
}
\arguments{
  \item{nrow}{number of rows in the grid, may also be an object inheriting from class \code{"SpatialGrid"} or \code{"GridTopology"} only in \code{cell2nb}}
  \item{ncol}{number of columns in the grid; if \code{nrow} or \code{x} is an object inheriting from class \code{"SpatialGrid"} or \code{"GridTopology"}, it may be omitted}
  \item{type}{default rook, may also be queen}
  \item{torus}{default FALSE, if TRUE map grid onto torus}
  \item{legacy}{default FALSE, nrow/ncol reversed, if TRUE wrong col/row directions (see \url{https://github.com/r-spatial/spdep/issues/20})}
  \item{x}{if given, an object inheriting from class \code{"SpatialGrid"} or \code{"GridTopology"}, and replaces \code{nrow} and \code{ncol}}
  \item{i}{vector of vector indices corresponding to rowcol, a matrix with two columns of row, column indices}
}

\value{
The function returns an object of class \code{nb} with a list of integer vectors containing neighbour region number ids. See \code{\link{card}} for details of \dQuote{nb} objects.
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{summary.nb}}, \code{\link{card}}}

\examples{
nb7rt <- cell2nb(7, 7)
summary(nb7rt)
xyc <- attr(nb7rt, "region.id")
xy <- matrix(as.integer(unlist(strsplit(xyc, ":"))), ncol=2, byrow=TRUE)
plot(nb7rt, xy)
nb7rt <- cell2nb(7, 7, torus=TRUE)
summary(nb7rt)
run <- FALSE
if (require("sp", quietly=TRUE)) run <- TRUE
if (run) {
# https://github.com/r-spatial/spdep/issues/20
GT <- GridTopology(c(1, 1), c(1, 1), c(10, 50))
SPix <- as(SpatialGrid(GT), "SpatialPixels")
nb_rook_cont <- poly2nb(as(SPix, "SpatialPolygons"), queen=FALSE)
nb_rook_dist <- dnearneigh(coordinates(SPix), 0, 1.01)
all.equal(nb_rook_cont, nb_rook_dist, check.attributes=FALSE)
## [1] TRUE
}
if (run) {
t.nb <- cell2nb(GT, type='rook', legacy=TRUE)
isTRUE(all.equal(nb_rook_cont, t.nb, check.attributes=FALSE))
## [1] FALSE
}
if (run) {
t.nb <- cell2nb(GT, type='rook')
isTRUE(all.equal(nb_rook_cont, t.nb, check.attributes=FALSE))
## [1] TRUE
}
if (run) {
# https://github.com/r-spatial/spdep/issues/55
# problem reported in issue caused by rep() cycling in unexpected order
GT <- GridTopology(c(1, 1), c(1, 1), c(22, 11))
SPix <- as(SpatialGrid(GT), "SpatialPixels")
nb_rook_cont <- poly2nb(as(SPix, "SpatialPolygons"), queen=FALSE)
nb_rook_dist <- dnearneigh(coordinates(SPix), 0, 1.01)
all.equal(nb_rook_cont, nb_rook_dist, check.attributes=FALSE)
}
if (run) {
t.nb <- cell2nb(GT, type='rook', legacy=TRUE)
isTRUE(all.equal(nb_rook_cont, t.nb, check.attributes=FALSE))
## [1] FALSE
}
if (run) {
t.nb <- cell2nb(GT, type='rook', legacy=FALSE)
isTRUE(all.equal(nb_rook_cont, t.nb, check.attributes=FALSE))
## [1] TRUE
}
}

\keyword{spatial}
