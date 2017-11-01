\name{similar.listw}
\alias{similar.listw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create symmetric similar weights lists}
\description{
  From Ord's 1975 paper, it is known that the Jacobian for SAR models may be found by "symmetrizing" by similarity (the eigenvalues of similar matrices are identical, so the Jacobian is too). This applies only to styles "W" and "S" with underlying symmetric binary neighbour relations or symmetric general neighbour relations (so no k-nearest neighbour relations). The function is invoked automatically within the SAR fitting functions, to call \code{eigen} on a symmetric matrix for the default eigen method, or to make it possible to use the Matrix method on weights that can be "symmetrized" in this way.
}
\usage{
similar.listw(listw)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
}

\value{
  a \code{listw} object
}
\references{Ord, J. K. 1975 Estimation methods for models of spatial interaction, \emph{Journal of the American Statistical Association}, 70, 120-126}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{lagsarlm}}, \code{\link{errorsarlm}}}
\examples{
data(oldcol)
COL.W <- nb2listw(COL.nb, style="W")
COL.S <- nb2listw(COL.nb, style="S")
sum(log(1 - 0.5 * eigenw(COL.W)))
sum(log(1 - 0.5 * eigenw(similar.listw(COL.W))))
W_J <- as(as_dsTMatrix_listw(similar.listw(COL.W)), "CsparseMatrix")
I <- as_dsCMatrix_I(dim(W_J)[1])
c(determinant(I - 0.5 * W_J, logarithm=TRUE)$modulus)
sum(log(1 - 0.5 * eigenw(COL.S)))
sum(log(1 - 0.5 * eigenw(similar.listw(COL.S))))
W_J <- as(as_dsTMatrix_listw(similar.listw(COL.S)), "CsparseMatrix")
c(determinant(I - 0.5 * W_J, logarithm=TRUE)$modulus)
}
\keyword{spatial}
