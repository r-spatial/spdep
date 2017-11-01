% Copyright 2002-8 by Roger S. Bivand
\name{globalG.test}
\alias{globalG.test}
\title{Global G test for spatial autocorrelation}
\description{
The global G statistic for spatial autocorrelation, complementing the local Gi LISA measures: \code{\link{localG}}.
}
\usage{
globalG.test(x, listw, zero.policy=NULL, alternative="greater",
 spChk=NULL, adjust.n=TRUE, B1correct=TRUE)
}
\arguments{
  \item{x}{a numeric vector the same length as the neighbours list in listw}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}; if a sequence of distance bands is to be used, it is recommended that the weights style be binary (one of \code{c("B", "C", "U")}).}
  \item{zero.policy}{default NULL, use global option value; if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA}
  \item{alternative}{a character string specifying the alternative hypothesis, must be one of "greater" (default), "less" or "two.sided".}
  \item{spChk}{should the data vector names be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use \code{get.spChkOption()}}
  \item{adjust.n}{default TRUE, if FALSE the number of observations is not adjusted for no-neighbour observations, if TRUE, the number of observations is adjusted}
  \item{B1correct}{default TRUE, if TRUE, the erratum referenced below: "On page 195, the coefficient of W2 in B1, (just below center of the page) should be 6,  not 3." is applied; if FALSE, 3 is used (as in CrimeStat IV)}
}

\value{
A list with class \code{htest} containing the following components:
  \item{statistic}{the value of the standard deviate of Moran's I.}
  \item{p.value}{the p-value of the test.}
  \item{estimate}{the value of the observed statistic, its expectation and variance.}
  \item{alternative}{a character string describing the alternative hypothesis.}
  \item{data.name}{a character string giving the name(s) of the data.}
}

\references{Getis. A, Ord, J. K. 1992 The analysis of spatial association by 
use of distance statistics, \emph{Geographical Analysis}, 24, p. 195; see
also Getis. A, Ord, J. K. 1993 Erratum, \emph{Geographical Analysis}, 25, 
p. 276.}
\author{Hisaji ONO \email{hi-ono@mn.xdsl.ne.jp} and Roger Bivand 
\email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{localG}}}

\examples{
if (require(rgdal, quietly=TRUE)) {
example(nc.sids, package="spData")
sidsrate79 <- (1000*nc.sids$SID79)/nc.sids$BIR79
dists <- c(10, 20, 30, 33, 40, 50, 60, 70, 80, 90, 100)
ndists <- length(dists)
ZG <- vector(mode="list", length=ndists)
names(ZG) <- as.character(dists)
milesxy <- cbind(nc.sids$east, nc.sids$north)
for (i in 1:ndists) {
  thisnb <- dnearneigh(milesxy, 0, dists[i])
  thislw <- nb2listw(thisnb, style="B", zero.policy=TRUE)
  ZG[[i]] <- globalG.test(sidsrate79, thislw, zero.policy=TRUE)
}
t(sapply(ZG, function(x) c(x$estimate[1], x$statistic, p.value=unname(x$p.value))))
for (i in 1:ndists) {
  thisnb <- dnearneigh(milesxy, 0, dists[i])
  thislw <- nb2listw(thisnb, style="B", zero.policy=TRUE)
  ZG[[i]] <- globalG.test(sidsrate79, thislw, zero.policy=TRUE, alternative="two.sided")
}
t(sapply(ZG, function(x) c(x$estimate[1], x$statistic, p.value=unname(x$p.value))))
}
}
\keyword{spatial}
