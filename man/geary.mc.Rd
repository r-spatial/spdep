% Copyright 2001-8 by Roger S. Bivand
\name{geary.mc}
\alias{geary.mc}
\title{Permutation test for Geary's C statistic}
\description{
 A permutation test for Geary's C statistic calculated by using nsim random permutations of x for the given spatial weighting scheme, to establish the rank of the observed statistic in relation to the nsim simulated values.
}
\usage{
geary.mc(x, listw, nsim, zero.policy=attr(listw, "zero.policy"), alternative="greater",
 spChk=NULL, adjust.n=TRUE, return_boot=FALSE, na.action=na.fail, scale=TRUE)
}
\arguments{
  \item{x}{a numeric vector the same length as the neighbours list in listw}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{nsim}{number of permutations}
  \item{zero.policy}{default \code{attr(listw, "zero.policy")} as set when \code{listw} was created, if attribute not set, use global option value; if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA}
  \item{alternative}{a character string specifying the alternative hypothesis, must be one of "greater" (default), or "less"; this reversal corresponds to that on \code{\link{geary.test}} described in the section on the output statistic value, based on Cliff and Ord 1973, p. 21 (changed 2011-04-11, thanks to Daniel Garavito).}
  \item{spChk}{should the data vector names be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use \code{get.spChkOption()}}
  \item{adjust.n}{default TRUE, if FALSE the number of observations is not adjusted for no-neighbour observations, if TRUE, the number of observations is adjusted}
  \item{return_boot}{return an object of class \code{boot} from the equivalent permutation bootstrap rather than an object of class \code{htest}}
  \item{na.action}{a function (default \code{na.fail}), can also be \code{na.omit} or \code{na.exclude} - in these cases the weights list will be subsetted to remove NAs in the data. It may be necessary to set zero.policy to TRUE because this subsetting may create no-neighbour observations. Note that only weights lists created without using the glist argument to \code{nb2listw} may be subsetted. \code{na.pass} is not permitted because it is meaningless in a permutation test.}
  \item{scale}{default TRUE, may be FALSE to revert changes made to accommodate \code{localC} in November 2021 (see #151)}
}

\value{
A list with class \code{htest} and \code{mc.sim} containing the following components:
  \item{statistic}{the value of the observed Geary's C.}
  \item{parameter}{the rank of the observed Geary's C.}
  \item{p.value}{the pseudo p-value of the test.}
  \item{alternative}{a character string describing the alternative hypothesis.}
  \item{method}{a character string giving the method used.}
  \item{data.name}{a character string giving the name(s) of the data, and the number of simulations.}
  \item{res}{nsim simulated values of statistic, final value is observed statistic}
}
\references{Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 63-5.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{geary}}, \code{\link{geary.test}}}

\examples{
data(oldcol)
set.seed(1)
sim1 <- geary.mc(COL.OLD$CRIME, nb2listw(COL.nb, style="W"),
 nsim=99, alternative="less")
sim1
mean(sim1$res)
var(sim1$res)
summary(sim1$res)
colold.lags <- nblag(COL.nb, 3)
sim2 <- geary.mc(COL.OLD$CRIME, nb2listw(colold.lags[[2]],
 style="W"), nsim=99)
sim2
summary(sim2$res)
sim3 <- geary.mc(COL.OLD$CRIME, nb2listw(colold.lags[[3]],
 style="W"), nsim=99)
sim3
summary(sim3$res)
crime <- COL.OLD$CRIME
is.na(crime) <- sample(1:length(crime), 10)
try(geary.mc(crime, nb2listw(COL.nb, style="W"), nsim=99,
 na.action=na.fail))
geary.mc(crime, nb2listw(COL.nb, style="W"), nsim=99, zero.policy=TRUE,
 na.action=na.omit)
geary.mc(crime, nb2listw(COL.nb, style="W"), nsim=99, zero.policy=TRUE,
 return_boot=TRUE, na.action=na.omit)
geary.mc(crime, nb2listw(COL.nb, style="W"), nsim=99, zero.policy=TRUE,
 na.action=na.exclude)
geary.mc(crime, nb2listw(COL.nb, style="W"), nsim=99, zero.policy=TRUE,
 return_boot=TRUE, na.action=na.exclude)
try(geary.mc(crime, nb2listw(COL.nb, style="W"), nsim=99, na.action=na.pass))
}
\keyword{spatial}
