% Copyright 2003-2010 by Roger S. Bivand
\name{mat2listw}
\alias{mat2listw}

\title{Convert a square spatial weights matrix to a weights list object}
\description{
  The function converts a square spatial weights matrix, 
optionally a sparse matrix to a weights list 
object, optionally adding region IDs from the row names of the matrix, as a 
sequence of numbers 1:nrow(x), or as given as an argument. The style can be imposed by rebuilting the weights list object internally.
}
\usage{
mat2listw(x, row.names = NULL, style="M")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A square non-negative matrix with no NAs representing spatial 
weights; may be a matrix of class \dQuote{sparseMatrix}}
  \item{row.names}{row names to use for region IDs}
  \item{style}{default "M", unknown style; if not "M", passed to \code{\link{nb2listw}} to re-build the object}
}

\value{
A \code{listw} object with the following members:
  \item{style}{"M", meaning matrix style, underlying style unknown, or assigned style argument in rebuilt object}
  \item{neighbours}{the derived neighbours list}
  \item{weights}{the weights for the neighbours derived from the matrix}
}

\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

 
\seealso{\code{\link{nb2listw}}, \code{\link{nb2mat}}}
\examples{
if (require(rgdal, quietly=TRUE)) {
example(columbus, package="spData")
coords <- coordinates(columbus)
col005 <- dnearneigh(coords, 0, 0.5, attr(col.gal.nb, "region.id"))
summary(col005)
col005.w.mat <- nb2mat(col005, zero.policy=TRUE)
col005.w.b <- mat2listw(col005.w.mat)
summary(col005.w.b$neighbours)
diffnb(col005, col005.w.b$neighbours)
col005.w.mat.3T <- kronecker(diag(3), col005.w.mat)
col005.w.b.3T <- mat2listw(col005.w.mat.3T, style="W")
summary(col005.w.b.3T$neighbours)
W <- as(nb2listw(col005, style="W", zero.policy=TRUE), "CsparseMatrix")
col005.spM <- mat2listw(W)
summary(col005.spM$neighbours)
diffnb(col005, col005.spM$neighbours)
IW <- kronecker(Diagonal(3), W)
col005.spM.3T <- mat2listw(IW, style="W")
summary(col005.spM.3T$neighbours)
}
}
\keyword{spatial}
