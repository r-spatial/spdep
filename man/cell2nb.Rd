% Copyright 2001 by Roger S. Bivand
\name{cell2nb}
\alias{cell2nb}
\alias{mrc2vi}
\alias{rookcell}
\alias{queencell}
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
cell2nb(nrow, ncol, type="rook", torus=FALSE)
mrc2vi(rowcol, nrow, ncol)
rookcell(rowcol, nrow, ncol, torus=FALSE, rmin=1, cmin=1)
queencell(rowcol, nrow, ncol, torus=FALSE, rmin=1, cmin=1)
vi2mrc(i, nrow, ncol)
}
\arguments{
  \item{nrow}{number of rows in the grid}
  \item{ncol}{number of columns in the grid}
  \item{type}{rook or queen}
  \item{torus}{map grid onto torus}
  \item{rowcol}{matrix with two columns of row, column indices}
  \item{i}{vector of vector indices corresponding to rowcol}
  \item{rmin}{lowest row index}
  \item{cmin}{lowset column index}
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
}
\keyword{spatial}
