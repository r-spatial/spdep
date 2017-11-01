% Copyright 2001 by Roger S. Bivand
\name{nb2mat}
\alias{nb2mat}
\alias{listw2mat}
\title{Spatial weights matrices for neighbours lists}
\description{
The function generates a weights matrix for a neighbours list with spatial weights for the chosen coding scheme.
}
\usage{
nb2mat(neighbours, glist=NULL, style="W", zero.policy=NULL)
listw2mat(listw)
}
\arguments{
  \item{neighbours}{an object of class \code{nb}}
  \item{glist}{list of general weights corresponding to neighbours}
  \item{style}{\code{style} can take values W, B, C, and S}
  \item{zero.policy}{default NULL, use global option value; if FALSE stop with error for any empty neighbour sets, if TRUE permit the weights list to be formed with zero-length weights vectors}
  \item{listw}{a \code{listw} object from for example \code{nb2listw}}
}
\details{
Starting from a binary neighbours list, in which regions are either listed as neighbours or are absent (thus not in the set of neighbours for some definition), the function creates an n by n weights matrix with values given by the coding scheme style chosen. B is the basic binary coding, W is row standardised, C is globally standardised, while S is the variance-stabilizing coding scheme proposed by Tiefelsdorf et al. 1999, p. 167-168.

The function leaves matrix rows as zero for any regions with zero neighbours fore zero.policy TRUE. These will in turn generate lag values of zero, equivalent to the sum of products of the zero row \code{t(rep(0, length=length(neighbours))) \%*\% x}, for arbitraty numerical vector \code{x} of length \code{length(neighbours)}. The spatially lagged value of x for the zero-neighbour region will then be zero, which may (or may not) be a sensible choice.
}
\value{
An n by n matrix, where n=length(neighbours)
}
\references{Tiefelsdorf, M., Griffith, D. A., Boots, B. 1999 A variance-stabilizing coding scheme for spatial link matrices, Environment and Planning A, 31, pp. 165-180.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{nb2listw}}}

\examples{
if (require(rgdal, quietly=TRUE)) {
example(columbus, package="spData")
coords <- coordinates(columbus)
col005 <- dnearneigh(coords, 0, 0.5, attr(col.gal.nb, "region.id"))
summary(col005)
col005.w.mat <- nb2mat(col005, zero.policy=TRUE)
table(round(apply(col005.w.mat, 1, sum)))
}}
\keyword{spatial}
