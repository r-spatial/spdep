\name{mstree}
\alias{mstree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find the minimal spanning tree}
\description{
  The minimal spanning tree is a connected graph with n nodes and n-1
  edges. This is a smaller class of possible partitions of a graph by
  pruning edges with high dissimilarity. If one edge is removed, the
  graph is partioned in two unconnected subgraphs. This function
  implements the algorithm due to Prim (1987).
}
\usage{
mstree(nbw, ini = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nbw}{An object of \code{listw} class returned by
    \code{\link{nb2listw}} function. See this help for details.}
  \item{ini}{The initial node in the minimal spanning tree.}
}
\details{
  The minimum spanning tree algorithm.

  Input a connected graph.

  Begin a empty set of nodes.

  Add an arbitrary note in this set.

  While are nodes not in the set, find a minimum cost edge connecting a
  node in the set and a node out of the set and add this node in the
  set.

  The set of edges is a minimum spanning tree.  
}
\value{
  A matrix with n-1 rows and tree columns. Each row is two nodes and the
  cost, i. e. the edge and it cost.
}
\references{
  R. C. Prim (1957) Shortest connection networks and some
  generalisations. In: Bell System Technical Journal, 36, pp. 1389-1401 
}
\author{Renato M. Assuncao and Elias T. Krainski}
%%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
### loading data
require(maptools)
bh <- readShapePoly(system.file("etc/shapes/bhicv.shp",
      package="spdep")[1])
### data padronized
dpad <- data.frame(scale(bh@data[,5:8]))

### neighboorhod list 
bh.nb <- poly2nb(bh)

### calculing costs
lcosts <- nbcosts(bh.nb, dpad)

### making listw
nb.w <- nb2listw(bh.nb, lcosts, style="B")

### find a minimum spanning tree
system.time(mst.bh <- mstree(nb.w,5))

dim(mst.bh)

head(mst.bh)
tail(mst.bh)

### the mstree plot
par(mar=c(0,0,0,0))
plot(mst.bh, coordinates(bh), col=2, 
     cex.lab=.7, cex.circles=0.035, fg="blue")
plot(bh, border=gray(.5), add=TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{graphs}
\keyword{spatial}% __ONLY ONE__ keyword per line
