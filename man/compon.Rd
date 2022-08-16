% Copyright 2001 by Roger S. Bivand and Nicholas Lewin-Koh
\name{Graph Components}
\alias{n.comp.nb}
%\alias{reach.ij}
%\alias{strong.comp}
\title{Depth First Search on Neighbor Lists}
\description{
  \code{n.comp.nb()} finds the number of disjoint connected subgraphs in the graph depicted by \code{nb.obj} - a spatial neighbours list object.
}
\usage{
n.comp.nb(nb.obj)
}
\arguments{
 \item{nb.obj}{a neighbours list object of class \code{nb}}
}
\value{
A list of:
\item{nc}{number of disjoint connected subgraphs}
\item{comp.id}{vector with the indices of the disjoint connected subgraphs that
 the nodes in \code{nb.obj} belong to}
}
\author{Nicholas Lewin-Koh \email{nikko@hailmail.net}}
\seealso{\code{\link{plot.nb}}}
\examples{
columbus <- st_read(system.file("shapes/columbus.shp", package="spData")[1], quiet=TRUE)
col.gal.nb <- read.gal(system.file("weights/columbus.gal", package="spData")[1])
coords <- st_coordinates(st_centroid(st_geometry(columbus)))
plot(col.gal.nb, coords, col="grey")
col2 <- droplinks(col.gal.nb, 21)
res <- n.comp.nb(col2)
table(res$comp.id)
plot(col2, coords, add=TRUE)
points(coords, col=res$comp.id, pch=16)
run <- FALSE
if (require(igraph, quietly=TRUE) && require(spatialreg, quietly=TRUE)) run <- TRUE
if (run) {
B <- as(nb2listw(col2, style="B", zero.policy=TRUE), "CsparseMatrix")
g1 <- graph.adjacency(B, mode="undirected")
c1 <- clusters(g1)
print(c1$no == res$nc)
}
if (run) {
print(all.equal(c1$membership, res$comp.id))
}
if (run) {
print(all.equal(c1$csize, c(table(res$comp.id)), check.attributes=FALSE))
}
if (run) {
W <- as(nb2listw(col2, style="W", zero.policy=TRUE), "CsparseMatrix")
g1W <- graph.adjacency(W, mode="directed", weighted="W")
c1W <- clusters(g1W)
print(all.equal(c1W$membership, res$comp.id, check.attributes=FALSE))
}
if (run) {
ow <- options("warn")$warn
options("warn"=2L)
# Matrix 1.4-2 vulnerability work-around
B1 <- try(get.adjacency(g1), silent=TRUE)
if (!inherits(B1, "try-error")) {
#B1 <- get.adjacency(g1)
print(all.equal(B, B1))
}
options("warn"=ow)
}
}

\keyword{spatial}
