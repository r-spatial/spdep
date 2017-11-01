\name{trW}
\Rdversion{1.1}
\alias{trW}
\alias{mom_calc}
\alias{mom_calc_int2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Spatial weights matrix powers traces}
\description{The function is used to prepare a vector of traces of powers of a spatial weights matrix}
\usage{
trW(W=NULL, m = 30, p = 16, type = "mult", listw=NULL, momentsSymmetry=TRUE)
mom_calc(lw, m)
mom_calc_int2(is, m, nb, weights, Card)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{W}{A spatial weights matrix in CsparseMatrix form}
  \item{m}{The number of powers; must be an even number for \sQuote{type}=\dQuote{moments} (default changed from 100 to 30 (2010-11-17))}
  \item{p}{The number of samples used in Monte Carlo simulation of the traces if type is MC (default changed from 50 to 16 (2010-11-17))}
  \item{type}{Either \dQuote{mult} (default) for powering a sparse matrix (with moderate or larger N, the matrix becomes dense, and may lead to swapping), or \dQuote{MC} for Monte Carlo simulation of the traces (the first two simulated traces are replaced by their analytical equivalents), or \dQuote{moments} to use the looping space saving algorithm proposed by Smirnov and Anselin (2009) - for \dQuote{moments}, \code{W} must be symmetric, for row-standardised weights through a similarity transformation}
  \item{listw, lw}{a listw object, which should either be fully symmetric, or be constructed as similar to symmetric from intrinsically symmetric neighbours using \code{\link{similar.listw}}, used with \sQuote{type}=\dQuote{moments}}
  \item{momentsSymmetry}{default TRUE; assert Smirnov/Anselin symmetry assumption}
  \item{is}{(used internally only in \code{mom_calc_int2} for \sQuote{type}=\dQuote{moments} on a cluster)}
  \item{nb}{(used internally only in \code{mom_calc_int2} for \sQuote{type}=\dQuote{moments} on a cluster)}
  \item{weights}{(used internally only in \code{mom_calc_int2} for \sQuote{type}=\dQuote{moments} on a cluster)}
  \item{Card}{(used internally only in \code{mom_calc_int2} for \sQuote{type}=\dQuote{moments} on a cluster)}
}

\value{
A numeric vector of \code{m} traces, with \dQuote{timings} and \dQuote{type} attributes; the \sQuote{type}=\dQuote{MC} also returns the standard deviation of the p-vector V divided by the square root of p as a measure of spread for the trace estimates.
}
\note{\code{mom_calc} and \code{mom_calc_int2} are for internal use only}
\references{LeSage J and RK Pace (2009) \emph{Introduction to Spatial Econometrics}. CRC Press, Boca Raton, pp. 96--105; Smirnov O and L Anselin (2009) An O(N) parallel method of computing the Log-Jacobian of the variable transformation for models with spatial interaction on a lattice. \emph{Computational Statistics and Data Analysis} 53 (2009) 2983--2984.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}


\seealso{
\code{\link{as_dgRMatrix_listw}}, \code{\link{nb2listw}}
}
\examples{
if (require(rgdal, quietly=TRUE)) {
example(columbus, package="spData")
listw <- nb2listw(col.gal.nb)
W <- as(listw, "CsparseMatrix")
system.time(trMat <- trW(W, type="mult"))
str(trMat)
set.seed(1100)
system.time(trMC <- trW(W, type="MC"))
str(trMC)
plot(trMat, trMC)
abline(a=0, b=1)
for(i in 3:length(trMC)) {
 segments(trMat[i], trMC[i]-2*attr(trMC, "sd")[i], trMat[i],
  trMC[i]+2*attr(trMC, "sd")[i])
}
listwS <- similar.listw(listw)
W <- forceSymmetric(as(listwS, "CsparseMatrix"))
system.time(trmom <- trW(W, m=24, type="moments"))
str(trmom)
all.equal(trMat[1:24], trmom, check.attributes=FALSE)
system.time(trMat <- trW(W, m=24, type="mult"))
str(trMat)
all.equal(trMat, trmom, check.attributes=FALSE)
set.seed(1)
system.time(trMC <- trW(W, m=24, type="MC"))
str(trMC)
}
\donttest{
data(boston, package="spData")
listw <- nb2listw(boston.soi)
listwS <- similar.listw(listw)
system.time(trmom <- trW(listw=listwS, m=24, type="moments"))
str(trmom)
library(parallel)
nc <- detectCores(logical=FALSE)
# set nc to 1L here
if (nc > 1L) nc <- 1L
coresOpt <- get.coresOption()
invisible(set.coresOption(nc))
if(!get.mcOption()) {
  cl <- makeCluster(get.coresOption())
  set.ClusterOption(cl)
}
system.time(trmomp <- trW(listw=listwS, m=24, type="moments"))
if(!get.mcOption()) {
  set.ClusterOption(NULL)
  stopCluster(cl)
}
all.equal(trmom, trmomp, check.attributes=FALSE)
invisible(set.coresOption(coresOpt))
}
}
\keyword{spatial}

