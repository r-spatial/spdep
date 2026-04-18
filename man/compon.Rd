% Copyright 2001-24 by Roger S. Bivand and Nicholas Lewin-Koh
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
\details{If \code{attr(nb.obj, "sym")} is \code{FALSE} and \code{igraph::components} is available, the components of the directed graph will be found by a simple breadth-first search; if \code{igraph::components} is not available, the object will be made symmetric (which may be time-consuming with large numbers of neighbours) and the components found by depth-first search. If \code{attr(nb.obj, "sym")} is \code{TRUE}, the components of the directed graph will be found by depth-first search. The time complexity of algorithms used in native code and through \code{igraph::components} is linear in the sum of the number of nodes and the number of edges in the graph, see \url{https://github.com/r-spatial/spdep/issues/160} for details; very dense neighbour objects will have large numbers of edges.}
\value{
A list of:
\item{nc}{number of disjoint connected subgraphs}
\item{comp.id}{vector with the indices of the disjoint connected subgraphs that
 the nodes in \code{nb.obj} belong to}
}
\author{Nicholas Lewin-Koh \email{nikko@hailmail.net}}
\seealso{\code{\link{plot.nb}}}
\examples{
columbus <- st_read(system.file("shapes/columbus.gpkg", package="spData")[1], quiet=TRUE)
col.gal.nb <- read.gal(system.file("weights/columbus.gal", package="spData")[1])
coords <- st_coordinates(st_centroid(st_geometry(columbus)))
plot(col.gal.nb, coords, col="grey")
col2 <- droplinks(col.gal.nb, 21)
res <- n.comp.nb(col2)
table(res$comp.id)
plot(col2, coords, add=TRUE)
points(coords, col=res$comp.id, pch=16)
run <- FALSE
if (require("igraph", quietly=TRUE) && require("spatialreg", quietly=TRUE)) run <- TRUE
if (run) {
B <- as(nb2listw(col2, style="B", zero.policy=TRUE), "CsparseMatrix")
g1 <- graph_from_adjacency_matrix(B, mode="undirected")
c1 <- components(g1)
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
g1W <- graph_from_adjacency_matrix(W, mode="directed", weighted="W")
c1W <- components(g1W, mode="weak")
print(all.equal(c1W$membership, res$comp.id, check.attributes=FALSE))
}
\donttest{
if (run) {
data(house, package="spData")
house <- sf::st_as_sf(house)
k6 <- knn2nb(knearneigh(house, k=6))
is.symmetric.nb(k6)
}
if (run) {
print(k6)
}
if (run) {
length(k6) + sum(card(k6))
}
if (run) {
# no pre-computed graph components
str(attr(k6, "ncomp"))
}
if (run) {
# raising the subgraph compute ceiling to above |N|+|E| computes and stores the
# object in the neighbour object
set.SubgraphCeiling(180000L)
k6 <- knn2nb(knearneigh(house, k=6))
str(attr(k6, "ncomp"))
}
if (run) {
print(k6)
}
if (run) {
system.time(udir <- n.comp.nb(make.sym.nb(k6)))
}
if (run) {
system.time(dir <- n.comp.nb(k6))
}
if (run) {
udir$nc
}
if (run) {
dir$nc
}
if (run) {
all.equal(dir, udir)
}
}
}


\keyword{spatial}
