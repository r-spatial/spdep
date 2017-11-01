% Copyright 2001 by Roger S. Bivand
\name{joincount.mc}
\alias{joincount.mc}

\title{Permutation test for same colour join count statistics}
\description{
A permutation test for same colour join count statistics calculated by using nsim random permutations of fx for the given spatial weighting scheme, to establish the ranks of the observed statistics (for each colour) in relation to the nsim simulated values.
}
\usage{
joincount.mc(fx, listw, nsim, zero.policy=FALSE, alternative="greater",
 spChk=NULL)
}
\arguments{
  \item{fx}{a factor of the same length as the neighbours and weights objects in listw}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{nsim}{number of permutations}
  \item{zero.policy}{if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA}
  \item{alternative}{a character string specifying the alternative hypothesis, must be one of "greater" (default), or "less".}
  \item{spChk}{should the data vector names be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use \code{get.spChkOption()}}
}

\value{
A list with class \code{jclist} of lists with class \code{htest} and \code{mc.sim} for each of the k colours containing the following components:
  \item{statistic}{the value of the observed statistic.}
  \item{parameter}{the rank of the observed statistic.}
  \item{method}{a character string giving the method used.}
  \item{data.name}{a character string giving the name(s) of the data.}
  \item{p.value}{the pseudo p-value of the test.}
  \item{alternative}{a character string describing the alternative hypothesis.}
  \item{estimate}{the mean and variance of the simulated distribution.}
  \item{res}{nsim simulated values of statistic, the final element is the observed statistic}
}
\references{Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 63-5.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{joincount.test}}}

\examples{
data(oldcol)
HICRIME <- cut(COL.OLD$CRIME, breaks=c(0,35,80), labels=c("low","high"))
names(HICRIME) <- rownames(COL.OLD)
joincount.mc(HICRIME, nb2listw(COL.nb, style="B"), nsim=99)
joincount.test(HICRIME, nb2listw(COL.nb, style="B"))
}
\keyword{spatial}
