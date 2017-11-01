\name{lextrB}
\alias{lextrB}
\alias{lextrW}
\alias{lextrS}
\alias{l_max}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find extreme eigenvalues of binary symmetric spatial weights}
\description{The functions find extreme eigenvalues of binary symmetric spatial weights, when these form planar graphs; general weights are not permiited. \code{l_max} finds the largest eigenvalue using Rayleigh quotient methods of any \dQuote{listw} object. \code{lextrB} first calls \code{l_max}, and uses its output to find the smallest eigenvalue in addition for binary symmetric spatial weights. \code{lextrW} extends these to find the smallest eigenvalue for intrinsically symmetric row-standardized binary weights matrices (transformed to symmetric through similarity internally). \code{lextrS} does the same for variance-stabilized (\dQuote{S} style) intrinsically symmetric binary weights matrices (transformed to symmetric through similarity internally).}
\usage{
lextrB(lw, zero.policy = TRUE, control = list())
lextrW(lw, zero.policy=TRUE, control=list())
lextrS(lw, zero.policy=TRUE, control=list())
l_max(lw, zero.policy=TRUE, control=list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lw}{a binary symmetric \code{listw} object from, for example, \code{nb2listw} with style \dQuote{B} for \code{lextrB}, style \dQuote{W} for \code{lextrW} and style \dQuote{S} for \code{lextrS}; for \code{l_max}, the object may be asymmetric and does not have to be binary}
  \item{zero.policy}{default NULL, use global option value; if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA}
  \item{control}{a list of control arguments}
}
\section{Control arguments}{
\describe{
  \item{trace}{report values in while loops, default NULL assuming FALSE; logical}
  \item{tol}{tolerance for breaking while loops, default \code{.Machine$double.eps^(1/2)}; numeric}
  \item{maxiter}{maximum number of iterations in while loops, default \code{6 * (length(lw$neighbours) - 2}; integer}
  \item{useC}{use C code, default TRUE, logical (not in \code{l_max})}
}
}
\value{
The functions return approximations to the extreme eigenvalues with the eigenvectors returned as attributes of this object.
}
\references{Griffith, D. A. (2004). Extreme eigenfunctions of adjacency matrices for planar graphs employed in spatial analyses. \emph{Linear Algebra and its Applications}, 388:201–219.}
\author{Roger Bivand, Yongwan Chun, Daniel Griffith}

\note{It may be necessary to modify control arguments if warnings about lack of convergence are seen.}

\examples{
data(boston, package="spData")
ab.listb <- nb2listw(boston.soi, style="B")
er <- range(eigenw(ab.listb))
er
res_1 <- lextrB(ab.listb)
c(res_1)
#if (require(igraph)) {
#  B <- as(ab.listb, "symmetricMatrix")
#  n <- length(boston.soi)
#  f2 <- function(x, extra=NULL) {as.vector(B \%*\% x)}
#  ar1 <- arpack(f2, sym=TRUE, options=list(n=n, nev=1, ncv=8,
#    which="LA", maxiter=200))
#  print(ar1$values)
# At line 409 of file dsaupd.f: Fortran runtime error: Actual string
# length is shorter than the declared one for dummy argument 'which' (0/2)
#  arn <- arpack(f2, sym=TRUE, options=list(n=n, nev=1, ncv=8,
#    which="SA", maxiter=200))
#  print(arn$values)
# At line 409 of file dsaupd.f: Fortran runtime error: Actual string
# length is shorter than the declared one for dummy argument 'which' (0/2)
#  ar1 <- arpack(f2, sym=TRUE, options=list(n=n, nev=2, ncv=8,
#    which="BE", maxiter=300))
# "BE" gives: At line 558 of file dsaup2.f: Fortran runtime error: 
# Index '9' of dimension 1 of array 'bounds' above upper bound of 8
# "BE" 
#  print(ar1$values)
#}
k5 <- knn2nb(knearneigh(boston.utm, k=5))
c(l_max(nb2listw(k5, style="B")))
max(Re(eigenw(nb2listw(k5, style="B"))))
c(l_max(nb2listw(k5, style="C")))
max(Re(eigenw(nb2listw(k5, style="C"))))
ab.listw <- nb2listw(boston.soi, style="W")
er <- range(eigenw(similar.listw(ab.listw)))
er
res_1 <- lextrW(ab.listw)
c(res_1)
#if (require(igraph)) {
#  B <- as(similar.listw(ab.listw), "symmetricMatrix")
#  ar1 <- arpack(f2, sym=TRUE, options=list(n=n, nev=1, ncv=8,
#    which="LA", maxiter=400))
#  print(ar1$values)
# At line 409 of file dsaupd.f: Fortran runtime error: Actual string
# length is shorter than the declared one for dummy argument 'which' (0/2)
#  arn <- arpack(f2, sym=TRUE, options=list(n=n, nev=1, ncv=8,
#    which="SA", maxiter=400))
#  print(arn$values)
# At line 409 of file dsaupd.f: Fortran runtime error: Actual string
# length is shorter than the declared one for dummy argument 'which' (0/2)
#  ar1 <- arpack(f2, sym=TRUE, options=list(n=n, nev=2, ncv=8,
#    which="BE", maxiter=300))
# "BE" gives: At line 558 of file dsaup2.f: Fortran runtime error: 
# Index '9' of dimension 1 of array 'bounds' above upper bound of 8
#  print(ar1$values)
#}
ab.listw <- nb2listw(boston.soi, style="S")
er <- range(eigenw(similar.listw(ab.listw)))
er
res_1 <- lextrS(ab.listw)
c(res_1)
#if (require(igraph)) {
#  B <- as(similar.listw(ab.listw), "symmetricMatrix")
#  ar1 <- arpack(f2, sym=TRUE, options=list(n=n, nev=1, ncv=8,
#    which="LA", maxiter=300))
#  print(ar1$values)
# At line 409 of file dsaupd.f: Fortran runtime error: Actual string
# length is shorter than the declared one for dummy argument 'which' (0/2)
#  arn <- arpack(f2, sym=TRUE, options=list(n=n, nev=1, ncv=8,
#    which="SA", maxiter=300))
#  print(arn$values)
# At line 409 of file dsaupd.f: Fortran runtime error: Actual string
# length is shorter than the declared one for dummy argument 'which' (0/2)
#  ar1 <- arpack(f2, sym=TRUE, options=list(n=n, nev=2, ncv=8,
#    which="BE", maxiter=300))
# "BE" gives: At line 558 of file dsaup2.f: Fortran runtime error: 
# Index '9' of dimension 1 of array 'bounds' above upper bound of 8
#  print(ar1$values)
#}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}
