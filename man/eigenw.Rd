% Copyright 2002-2012 by Roger S. Bivand
\name{eigenw}
\alias{eigenw}
\alias{griffith_sone}
\alias{subgraph_eigenw}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Spatial weights matrix eigenvalues}
\description{
  The \code{eigenw} function returns a numeric vector of eigenvalues of 
the weights matrix generated from the spatial weights object \code{listw}. 
The eigenvalues are used to speed the computation of the Jacobian in 
spatial model estimation:

\deqn{\log(\det[I - \rho W]) = \sum_{i=1}^{n}\log(1 - \rho \lambda_i)}

where \eqn{W}{W} is the n by n spatial weights matrix, and \eqn{\lambda_i}{lambda[i]} are the
eigenvalues of \eqn{W}{W}.
}
\usage{
eigenw(listw, quiet=NULL)
griffith_sone(P, Q, type="rook")
subgraph_eigenw(nb, glist=NULL, style="W", zero.policy=NULL, quiet=NULL)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{quiet}{default NULL, use global !verbose option value; set to FALSE for short summary}
  \item{P}{number of columns in the grid (number of units in a horizontal axis direction)}
  \item{Q}{number of rows in the grid (number of units in a vertical axis direction.)}
  \item{type}{\dQuote{rook} or \dQuote{queen}}
  \item{nb}{an object of class \code{nb}}
  \item{glist}{list of general weights corresponding to neighbours}
  \item{style}{\code{style} can take values \dQuote{W}, \dQuote{B}, \dQuote{C}, \dQuote{U}, \dQuote{minmax} and \dQuote{S}}
  \item{zero.policy}{default NULL, use global option value; if FALSE stop with error for any empty neighbour sets, if TRUE permit the weights list to be formed with zero-length weights vectors}
}

\details{The \code{eigenw} function computes the eigenvalues of a single spatial weights object. The \code{griffith_sone} function function may be used, following Ord and Gasim (for references see Griffith and Sone (1995)), to calculate analytical eigenvalues for binary rook or queen contiguous neighbours where the data are arranged as a regular P times Q grid. The \code{subgraph_eigenw} function may be used when there are multiple graph components, of which the largest may be handled as a dense matrix. Here the eigenvalues are computed for each subgraph in turn, and catenated to reconstruct the complete set. The functions may be used to provide pre-computed eigenvalues for spatial regression functions.}

\value{
  a numeric or complex vector of eigenvalues of the weights matrix generated from the spatial weights object.
}
\references{Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 155;
Ord, J. K. 1975 Estimation methods for models of spatial interaction, Journal
of the American Statistical Association, 70, 120-126.; Griffith, D. A. and Sone, A. (1995). Trade-offs associated with normalizing constant computational simplifications for estimating spatial statistical models. Journal of Statistical Computation and Simulation, 51, 165-183.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{eigen}}%, \code{\link{logSpwdet}}
}

\examples{
data(oldcol)
W.eig <- eigenw(nb2listw(COL.nb, style="W"))
1/range(W.eig)
S.eig <- eigenw(nb2listw(COL.nb, style="S"))
1/range(S.eig)
B.eig <- eigenw(nb2listw(COL.nb, style="B"))
1/range(B.eig)
# cases for intrinsically asymmetric weights
crds <- cbind(COL.OLD$X, COL.OLD$Y)
k3 <- knn2nb(knearneigh(crds, k=3))
is.symmetric.nb(k3)
k3eig <- eigenw(nb2listw(k3, style="W"))
is.complex(k3eig)
rho <- 0.5
Jc <- sum(log(1 - rho * k3eig))
# complex eigenvalue Jacobian
Jc
# subgraphs
nc <- n.comp.nb(k3)
nc$nc
table(nc$comp.id)
k3eigSG <- subgraph_eigenw(k3, style="W")
all.equal(sort(k3eig), k3eigSG)
W <- as(nb2listw(k3, style="W"), "CsparseMatrix")
I <- diag(length(k3))
Jl <- sum(log(abs(diag(slot(lu(I - rho * W), "U")))))
# LU Jacobian equals complex eigenvalue Jacobian
Jl
all.equal(Re(Jc), Jl)
# wrong value if only real part used
Jr <- sum(log(1 - rho * Re(k3eig)))
Jr
all.equal(Jr, Jl)
# construction of Jacobian from complex conjugate pairs (Jan Hauke)
Rev <- Re(k3eig)[which(Im(k3eig) == 0)]
# real eigenvalues
Cev <- k3eig[which(Im(k3eig) != 0)]
pCev <- Cev[Im(Cev) > 0]
# separate complex conjugate pairs
RpCev <- Re(pCev)
IpCev <- Im(pCev)
# reassemble Jacobian
Jc1 <- sum(log(1 - rho*Rev)) + sum(log((1 - rho * RpCev)^2 + (rho^2)*(IpCev^2)))
all.equal(Re(Jc), Jc1)
# impact of omitted complex part term in real part only Jacobian
Jc2 <- sum(log(1 - rho*Rev)) + sum(log((1 - rho * RpCev)^2))
all.equal(Jr, Jc2)
# trace of asymmetric (WW) and crossprod of complex eigenvalues for APLE
sum(diag(W \%*\% W))
crossprod(k3eig)
# analytical regular grid eigenvalues
rg <- cell2nb(ncol=7, nrow=7, type="rook")
rg_eig <- eigenw(nb2listw(rg, style="B"))
rg_GS <- griffith_sone(P=7, Q=7, type="rook")
all.equal(rg_eig, rg_GS)
\donttest{
if (require(igraph)) {
B <- as(nb2listw(rg, style="B"), "CsparseMatrix")
f <- function(x, extra=NULL) {as.vector(B \%*\% x)}
res1 <- arpack(f, sym=TRUE, options=list(n=nrow(B), nev=1, ncv=8, which="LA",
 maxiter=200))
resn <- arpack(f, sym=TRUE, options=list(n=nrow(B), nev=1, ncv=8, which="SA", maxiter=200))
print(c(resn$value, res1$value))
#res <- arpack(f, sym=TRUE, options=list(n=nrow(B), nev=2, ncv=8, which="BE",
# maxiter=200))
#print(res$value)
# At line 558 of file dsaup2.f: Fortran runtime error:
# Index '9' of dimension 1 of array 'bounds' above upper bound of 8
print(all.equal(range(Re(rg_eig)), c(resn$value, res1$value)))
lw <- nb2listw(rg, style="W")
rg_eig <- eigenw(similar.listw(lw))
print(range(Re(rg_eig)))
W  <- as(lw, "CsparseMatrix")
f <- function(x, extra=NULL) {as.vector(W \%*\% x)}
print(arpack(f, sym=FALSE, options=list(n=nrow(W), nev=1, ncv=8, which="LR",
 maxiter=200))$value)
print(arpack(f, sym=FALSE, options=list(n=nrow(W), nev=1, ncv=8, which="SR",
 maxiter=200))$value)
}}
}
\keyword{spatial}
