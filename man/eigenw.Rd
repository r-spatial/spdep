% Copyright 2002-2012 by Roger S. Bivand
\name{eigenw}
\alias{eigenw}
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
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{quiet}{default NULL, use global !verbose option value; set to FALSE for short summary}
}

\details{The \code{eigenw} function computes the eigenvalues of a single spatial weights object. }

\value{
  a numeric or complex vector of eigenvalues of the weights matrix generated from the spatial weights object.
}
\references{Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 155;
Ord, J. K. 1975 Estimation methods for models of spatial interaction, Journal
of the American Statistical Association, 70, 120-126.}
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
\dontrun{
run <- FALSE
if (require("RSpectra", quietly=TRUE)) run <- TRUE
if (run) {
rg <- cell2nb(ncol=7, nrow=7, type="rook")
B <- as(nb2listw(rg, style="B"), "CsparseMatrix")
res1 <- eigs(B, k=1, which="LR")$values
resn <- eigs(B, k=1, which="SR")$values
print(Re(c(resn, res1)))
}
if (run) {
rg_eig <- eigenw(nb2listw(rg, style="B"))
print(all.equal(range(Re(rg_eig)), c(resn, res1))) 
}
if (run) {
lw <- nb2listw(rg, style="W")
rg_eig <- eigenw(similar.listw(lw))
print(range(Re(rg_eig)))
}
if (run) {
W  <- as(lw, "CsparseMatrix")
print(Re(c(eigs(W, k=1, which="SR")$values, eigs(W, k=1, which="LR")$values)))
}}
}
\keyword{spatial}
