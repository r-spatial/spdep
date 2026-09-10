\name{normalized.moran.mc}
\alias{normalized.moran.mc}
\title{Permutation test for normalized Moran's \eqn{I} statistic}
\usage{
normalized.moran.mc(
  x,
  listw,
  nsim,
  zero.policy = attr(listw, "zero.policy"),
  alternative = "greater",
  na.action = na.fail,
  spChk = NULL,
  return_boot = FALSE,
  adjust.n = TRUE
)
}
\arguments{
\item{x}{a numeric vector the same length as the neighbours list in \code{listw}}

\item{listw}{a \code{listw} object created for example by nb2listw}

\item{nsim}{number of permutations}

\item{zero.policy}{default \code{attr(listw, "zero.policy")} as set when \code{listw} was created, if attribute not set, use global option value; if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA}

\item{alternative}{a character string specifying the alternative hypothesis, must be one of \code{greater} (default), \code{less} or \code{two.sided}}

\item{na.action}{a function (default \code{na.fail}), can also be \code{na.omit} or \code{na.exclude} - in these cases the weights list will be subsetted to remove NAs in the data. It may be necessary to set \code{zero.policy} to TRUE because this subsetting may create no-neighbour observations. Note that only weights lists created without using the \code{glist} argument to \code{nb2listw} may be subsetted. If \code{na.pass} is used, zero is substituted for NA values in calculating the spatial lag}

\item{spChk}{should the data vector names be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use \code{get.spChkOption()}}

\item{return_boot}{return an object of class \code{boot} from the equivalent permutation bootstrap rather than an object of class \code{htest}}

\item{adjust.n}{default TRUE, if FALSE the number of observations is not adjusted for no-neighbour observations, if TRUE, the number of observations is adjusted}
}
\value{
A list with class \code{htest} and \code{mc.sim} containing the following components: \describe{
\item{\code{statistic}}{the value of the observed normalized Moran's I.}
\item{\code{parameter}}{the rank of the observed normalized Moran's I.}
\item{\code{p.value}}{the pseudo p-value of the test.}
\item{\code{alternative}}{a character string describing the alternative hypothesis.}
\item{\code{method}}{a character string giving the method used.}
\item{\code{data.name}}{a character string giving the name(s) of the data.}
\item{\code{res}}{\code{nsim} simulated values of statistic, final value is observed statistic.}
}
}
\description{
A permutation test for the "\eqn{I_6}" normalized Moran's \eqn{I} statistic by Tillé et al. (2026) calculated by using \code{nsim} random permutations of \code{x} for the given spatial weighting scheme, to establish the rank of the observed statistic in relation to the \code{nsim} simulated values.
}
\section{Author}{
 Diego Giuliani \email{diego.giuliani@unitn.it}
}

\section{References}{
 Tillé Y., Giuliani D., Dickson M.M., Espa G. (2026) Revisiting the Normalization of the Moran's \eqn{I} index: A Correlation-based Approach with Inference, \emph{International Statistical Review}
}

\examples{
library(spdep)
data(oldcol)
colw <- nb2listw(COL.nb, style="W")
nsim <- 99
set.seed(1234)
sim1 <- normalized.moran.mc(COL.OLD$CRIME, listw=colw, nsim=nsim)
sim1
mean(sim1$res[1:nsim])
var(sim1$res[1:nsim])
summary(sim1$res[1:nsim])
colold.lags <- nblag(COL.nb, 3)
set.seed(1234)
sim2 <- normalized.moran.mc(COL.OLD$CRIME, nb2listw(colold.lags[[2]],
                                                    style="W"), nsim=nsim)
summary(sim2$res[1:nsim])
sim3 <- normalized.moran.mc(COL.OLD$CRIME, nb2listw(colold.lags[[3]],
                                                    style="W"), nsim=nsim)
summary(sim3$res[1:nsim])
crime <- COL.OLD$CRIME
is.na(crime) <- sample(1:length(crime), 10)
try(normalized.moran.mc(crime, nb2listw(COL.nb, style="W"), nsim=99,
                        na.action=na.fail))
normalized.moran.mc(crime, nb2listw(COL.nb, style="W"), nsim=99, zero.policy=TRUE,
                    na.action=na.omit)
normalized.moran.mc(crime, nb2listw(COL.nb, style="W"), nsim=99, zero.policy=TRUE,
                    return_boot=TRUE, na.action=na.omit)
normalized.moran.mc(crime, nb2listw(COL.nb, style="W"), nsim=99, zero.policy=TRUE,
                    na.action=na.exclude)
normalized.moran.mc(crime, nb2listw(COL.nb, style="W"), nsim=99, zero.policy=TRUE,
                    return_boot=TRUE, na.action=na.exclude)
try(normalized.moran.mc(crime, nb2listw(COL.nb, style="W"), nsim=99, na.action=na.pass))
}
\seealso{
\link{normalized.moran}, \link{normalized.moran.test}
}
