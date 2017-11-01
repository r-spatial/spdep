% Copyright 2001-9 by Roger S. Bivand
\name{invIrM}
\alias{invIrM}
\alias{invIrW}
\alias{powerWeights}

\title{Compute SAR generating operator}
\description{
  Computes the matrix used for generating simultaneous autoregressive random variables, for a given value of rho, a neighbours list object or a matrix, a chosen coding scheme style, and optionally a list of general weights corresponding to neighbours.
}
\usage{
invIrM(neighbours, rho, glist=NULL, style="W", method="solve",
 feasible=NULL)
invIrW(x, rho, method="solve", feasible=NULL)
powerWeights(W, rho, order=250, X, tol=.Machine$double.eps^(3/5))
}

\arguments{
  \item{neighbours}{an object of class \code{nb}}
  \item{rho}{autoregressive parameter}
  \item{glist}{list of general weights corresponding to neighbours}
  \item{style}{\code{style} can take values W, B, C, and S}
  \item{method}{default \code{solve}, can also take value \code{chol}}
  \item{feasible}{if NULL, the given value of rho is checked to see if it lies within its feasible range, if TRUE, the test is not conducted}
  \item{x}{either a \code{listw} object from for example \code{nb2listw}, or a square spatial weights matrix, optionally a sparse matrix}
  \item{W}{A spatial weights matrix (either a dense matrix or a CsparseMatrix)}
  \item{order}{Power series maximum limit}
  \item{X}{A numerical matrix}
  \item{tol}{Tolerance for convergence of power series}
}
\details{
The \code{invIrW} function generates the full weights matrix V, checks that rho lies in its feasible range between 1/min(eigen(V)) and 1/max(eigen(V)), and returns the nxn inverted matrix \deqn{(I - \rho V)^{-1}}. With method=\dQuote{chol} (only for a listw object), Cholesky decomposition is used, thanks to contributed code by Markus Reder and Werner Mueller.

The \code{powerWeights} function uses power series summation to cumulate the product \deqn{(I - \rho V)^{-1} \%*\% X} from \deqn{(I + \rho V + (\rho V)^2 + \dots) \%*\% X}, which can be done by storing only sparse V and several matrices of the same dimensions as X. This makes it possible to handle larger spatial weights matrices, but is sensitive to the power weights order and the tolerance arguments when the spatial coefficient is close to its bounds, leading to incorrect estimates of the implied inverse matrix.
}
\value{
  An nxn matrix with a "call" attribute; the \code{powerWeights} function returns a matrix of the same dimensions as X which has been multipled by the power series equivalent of the dense matrix \deqn{(I - \rho V)^{-1}}.
}
\note{
Before version 0.6-10, \code{powerWeights} only worked correctly for positive rho, with differences from true values increasing as rho approached -1, and exploding between -1 and the true negative bound.
}
\references{Tiefelsdorf, M., Griffith, D. A., Boots, B. 1999 A variance-stabilizing coding scheme for spatial link matrices, Environment and Planning A, 31, pp. 165-180; Tiefelsdorf, M. 2000 Modelling spatial processes, Lecture notes in earth sciences, Springer, p. 76; Haining, R. 1990 Spatial data analysis in the social and environmental sciences, Cambridge University Press, p. 117; Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 152; Reder, M. and Mueller, W. (2007) An Improvement of the invIrM Routine of the Geostatistical R-package spdep by Cholesky Inversion, Statistical Projects, LV No: 238.205, SS 2006, Department of Applied Statistics, Johannes Kepler University, Linz}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{nb2listw}}}

\examples{

nb7rt <- cell2nb(7, 7, torus=TRUE)
set.seed(1)
x <- matrix(rnorm(500*length(nb7rt)), nrow=length(nb7rt))
res0 <- apply(invIrM(nb7rt, rho=0.0, method="chol",
 feasible=TRUE) \%*\% x, 2, function(x) var(x)/length(x))
res2 <- apply(invIrM(nb7rt, rho=0.2, method="chol",
 feasible=TRUE) \%*\% x, 2, function(x) var(x)/length(x))
res4 <- apply(invIrM(nb7rt, rho=0.4, method="chol",
 feasible=TRUE) \%*\% x, 2, function(x) var(x)/length(x))
res6 <- apply(invIrM(nb7rt, rho=0.6, method="chol",
 feasible=TRUE) \%*\% x, 2, function(x) var(x)/length(x))
res8 <- apply(invIrM(nb7rt, rho=0.8, method="chol",
 feasible=TRUE) \%*\% x, 2, function(x) var(x)/length(x))
res9 <- apply(invIrM(nb7rt, rho=0.9, method="chol",
 feasible=TRUE) \%*\% x, 2, function(x) var(x)/length(x))
plot(density(res9), col="red", xlim=c(-0.01, max(density(res9)$x)),
  ylim=range(density(res0)$y),
  xlab="estimated variance of the mean",
  main=expression(paste("Effects of spatial autocorrelation for different ",
    rho, " values")))
lines(density(res0), col="black")
lines(density(res2), col="brown")
lines(density(res4), col="green")
lines(density(res6), col="orange")
lines(density(res8), col="pink")
legend(c(-0.02, 0.01), c(7, 25),
 legend=c("0.0", "0.2", "0.4", "0.6", "0.8", "0.9"),
 col=c("black", "brown", "green", "orange", "pink", "red"), lty=1, bty="n")
\dontrun{
x <- matrix(rnorm(length(nb7rt)), ncol=1)
system.time(e <- invIrM(nb7rt, rho=0.9, method="chol", feasible=TRUE) \%*\% x)
system.time(e <- invIrM(nb7rt, rho=0.9, method="chol", feasible=NULL) \%*\% x)
system.time(e <- invIrM(nb7rt, rho=0.9, method="solve", feasible=TRUE) \%*\% x)
system.time(e <- invIrM(nb7rt, rho=0.9, method="solve", feasible=NULL) \%*\% x)
W <- as(nb2listw(nb7rt), "CsparseMatrix")
system.time(ee <- powerWeights(W, rho=0.9, X=x))
all.equal(e, as(ee, "matrix"), check.attributes=FALSE)
system.time(e <- invIrM(nb7rt, rho=0.98, method="solve", feasible=NULL) \%*\% x)
system.time(ee <- powerWeights(W, rho=0.98, X=x))
str(attr(ee, "internal"))
all.equal(e, as(ee, "matrix"), check.attributes=FALSE)
system.time(ee <- powerWeights(W, rho=0.98, order=1000, X=x))
all.equal(e, as(ee, "matrix"), check.attributes=FALSE)
nb60rt <- cell2nb(60, 60, torus=TRUE)
W <- as(nb2listw(nb60rt), "CsparseMatrix")
set.seed(1)
x <- matrix(rnorm(dim(W)[1]), ncol=1)
system.time(ee <- powerWeights(W, rho=0.3, X=x))
str(as(ee, "matrix"))
obj <- errorsarlm(as(ee, "matrix")[,1] ~ 1, listw=nb2listw(nb60rt), method="Matrix")
coefficients(obj)
}
}
\keyword{spatial}
