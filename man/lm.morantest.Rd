% Copyright 2001-5 by Roger S. Bivand
\name{lm.morantest}
\alias{lm.morantest}
\alias{listw2U}
\title{Moran's I test for residual spatial autocorrelation}
\description{
  Moran's I test for spatial autocorrelation in residuals from an estimated linear model (\code{lm()}). The helper function \code{listw2U()} constructs a weights list object corresponding to the sparse matrix \eqn{\frac{1}{2} ( \mathbf{W} + \mathbf{W}'}{1/2 (W + W')}
}
\usage{
lm.morantest(model, listw, zero.policy=NULL, alternative = "greater",
 spChk=NULL, resfun=weighted.residuals, naSubset=TRUE)
listw2U(listw)
}
\arguments{
  \item{model}{an object of class \code{lm} returned by \code{lm}; weights
may be specified in the \code{lm} fit, but offsets should not be used}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{zero.policy}{default NULL, use global option value; if TRUE assign zero to the lagged value of zones without
neighbours, if FALSE assign NA}
  \item{alternative}{a character string specifying the alternative hypothesis,
must be one of "greater" (default), "less" or "two.sided".}
  \item{spChk}{should the data vector names be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use \code{get.spChkOption()}}
  \item{resfun}{default: weighted.residuals; the function to be used to extract residuals from the \code{lm} object, may be \code{residuals}, \code{weighted.residuals}, \code{rstandard}, or \code{rstudent}}
  \item{naSubset}{default TRUE to subset listw object for omitted observations in model object (this is a change from earlier behaviour, when the \code{model$na.action} component was ignored, and the listw object had to be subsetted by hand)}
}
\value{
A list with class \code{htest} containing the following components:
  \item{statistic}{the value of the standard deviate of Moran's I.}
  \item{p.value}{the p-value of the test.}
  \item{estimate}{the value of the observed Moran's I, its expectation and
variance under the method assumption.}
  \item{alternative}{a character string describing the alternative hypothesis.}
  \item{method}{a character string giving the method used.}
  \item{data.name}{a character string giving the name(s) of the data.}
}
\references{Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 203, }
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{lm.LMtests}}, \code{\link{lm}}}

\examples{
data(oldcol)
oldcrime1.lm <- lm(CRIME ~ 1, data = COL.OLD)
oldcrime.lm <- lm(CRIME ~ HOVAL + INC, data = COL.OLD)
lm.morantest(oldcrime.lm, nb2listw(COL.nb, style="W"))
lm.LMtests(oldcrime.lm, nb2listw(COL.nb, style="W"))
lm.morantest(oldcrime.lm, nb2listw(COL.nb, style="S"))
lm.morantest(oldcrime1.lm, nb2listw(COL.nb, style="W"))
moran.test(COL.OLD$CRIME, nb2listw(COL.nb, style="W"),
 randomisation=FALSE)
oldcrime.wlm <- lm(CRIME ~ HOVAL + INC, data = COL.OLD,
 weights = I(1/AREA_PL))
lm.morantest(oldcrime.wlm, nb2listw(COL.nb, style="W"),
 resfun=weighted.residuals)
lm.morantest(oldcrime.wlm, nb2listw(COL.nb, style="W"),
 resfun=rstudent)
}
\keyword{spatial}
