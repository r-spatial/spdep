% Copyright 2002-8 by Roger S. Bivand
\name{globalG.test}
\alias{globalG.test}
\title{Global G test for spatial autocorrelation}
\description{
The global G statistic for spatial autocorrelation, complementing the local Gi LISA measures: \code{\link{localG}}.
}
\usage{
globalG.test(x, listw, zero.policy=attr(listw, "zero.policy"), alternative="greater",
 spChk=NULL, adjust.n=TRUE, B1correct=TRUE, adjust.x=TRUE, Arc_all_x=FALSE,
 na.action=na.fail)
}
\arguments{
  \item{x}{a numeric vector the same length as the neighbours list in listw}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}; if a sequence of distance bands is to be used, it is recommended that the weights style be binary (one of \code{c("B", "C", "U")}).}
  \item{zero.policy}{default \code{attr(listw, "zero.policy")} as set when \code{listw} was created, if attribute not set, use global option value; if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA}
  \item{alternative}{a character string specifying the alternative hypothesis, must be one of "greater" (default), "less" or "two.sided".}
  \item{spChk}{should the data vector names be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use \code{get.spChkOption()}}
  \item{adjust.n}{default TRUE, if FALSE the number of observations is not adjusted for no-neighbour observations, if TRUE, the number of observations is adjusted}
  \item{B1correct}{default TRUE, if TRUE, the erratum referenced below: "On page 195, the coefficient of W2 in B1, (just below center of the page) should be 6,  not 3." is applied; if FALSE, 3 is used (as in CrimeStat IV)}
  \item{adjust.x}{default TRUE, if TRUE, x values of observations with no neighbours are omitted in the denominator of G}
  \item{Arc_all_x}{default FALSE, if Arc_all_x=TRUE and adjust.x=TRUE, use the full x vector in part of the denominator term for G}
  \item{na.action}{a function (default \code{na.fail}), can also be \code{na.omit} or \code{na.exclude} - in these cases the weights list will be subsetted to remove NAs in the data. It may be necessary to set zero.policy to TRUE because this subsetting may create no-neighbour observations. Note that only weights lists created without using the glist argument to \code{nb2listw} may be subsetted. \code{na.pass} is not permitted.}
}

\value{
A list with class \code{htest} containing the following components:
  \item{statistic}{the value of the standard deviate of Getis-Ord G.}
  \item{p.value}{the p-value of the test.}
  \item{estimate}{the value of the observed statistic, its expectation and variance.}
  \item{alternative}{a character string describing the alternative hypothesis.}
  \item{data.name}{a character string giving the name(s) of the data.}
}

\references{Getis. A, Ord, J. K. 1992 The analysis of spatial association by 
use of distance statistics, \emph{Geographical Analysis}, 24, p. 195; see
also Getis. A, Ord, J. K. 1993 Erratum, \emph{Geographical Analysis}, 25, 
p. 276; Bivand RS, Wong DWS 2018 Comparing implementations of global and local indicators of spatial association. TEST, 27(3), 716--748 \doi{10.1007/s11749-018-0599-x}}
\author{Hisaji ONO \email{hi-ono@mn.xdsl.ne.jp} and Roger Bivand 
\email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{localG}}}

\examples{
nc.sids <- st_read(system.file("shapes/sids.gpkg", package="spData")[1], quiet=TRUE)
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
data(oldcol)
crime <- COL.OLD$CRIME
is.na(crime) <- sample(1:length(crime), 10)
res <- try(globalG.test(crime, nb2listw(COL.nb, style="B"),
 na.action=na.fail))
globalG.test(crime, nb2listw(COL.nb, style="B"), zero.policy=TRUE,
 na.action=na.omit)
globalG.test(crime, nb2listw(COL.nb, style="B"), zero.policy=TRUE,
 na.action=na.exclude)
try(globalG.test(crime, nb2listw(COL.nb, style="B"), na.action=na.pass))
}
\keyword{spatial}
