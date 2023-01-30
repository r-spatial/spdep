\name{grid2nb}
\alias{grid2nb}
\title{Construct neighbours for a GridTopology}
\description{
  The function builds a neighbours list for a grid topology.
  It works for a k-dimentional grid topology, k>=1.}
\usage{
grid2nb(grid, d = grid@cells.dim,
        queen = TRUE, nb = TRUE, self = FALSE)
}

\arguments{
  \item{grid}{An object of class \code{GridTopology}. One can avoid to 
  supply this by just suplying the dimentions in the \code{d} argument.} 
  \item{d}{A scalar (for one dimentional grid) or a length k vector
  specyfying the number of grid cells in each direction of the k
  dimentions.}
  \item{queen}{Logical. Default is TRUE. To inform if the queen
  neighbourhood structure should be considered. If FALSE, only 
  a hyper-cube with a common face will be considered neighbour. 
  If TRUE, a single shared coordinate meets the contiguity condition.} 
  \item{nb}{Default TRUE. If TRUE, return the result as a neighbours
  list with class \code{nb}. If FALSE, the result is a matrix with
  \code{3^k} columns if \code{self = TRUE} or \code{3^k-1} if
  \code{self = FALSE}. Zeros are used for hyper-cubes at boundaries.}
  \item{self}{Default FALSE, to indicate if the hyper-cube neighbour
  itself should be considered a neighbour.}
}
\value{
  Either a matrix, if \dQuote{nb} is FALSE or a neighbours list with
  class \code{nb}. See \code{\link{card}} for details of \dQuote{nb}
  objects.}

\note{
  This applies to a k-dimentional grid topology.
}
\author{Elias T Krainski \email{eliaskrainski@gmail.com}}

\seealso{\code{\link{poly2nb}}, \code{\link{summary.nb}},
  \code{\link{card}}} 

\examples{
nb <- grid2nb(d = c(5L, 5L, 5L))
nb
summary(nb)
if (require("sp", quietly=TRUE)) {
gt <- GridTopology(c(.125,.1), c(.25,.2), c(4L, 5L))
nb1 <- grid2nb(gt, queen = FALSE)
nb2 <- grid2nb(gt)

sg <- SpatialGrid(gt)
plot(sg, lwd=3)
plot(nb1, coordinates(sg), add=TRUE, lty=2, col=2, lwd=2)
plot(nb2, coordinates(sg), add=TRUE, lty=3, col=4)

str(grid2nb(d=5))
}
}

\keyword{spatial}
