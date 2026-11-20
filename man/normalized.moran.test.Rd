\name{normalized.moran.test}
\alias{normalized.moran.test}
\title{Normalized Moran's \eqn{I} test for spatial autocorrelation}
\usage{
normalized.moran.test(
  x,
  listw,
  zero.policy = attr(listw, "zero.policy"),
  alternative = "greater",
  na.action = na.fail,
  spChk = NULL,
  adjust.n = TRUE
)
}
\arguments{
\item{x}{a numeric vector the same length as the neighbours list in listw}

\item{listw}{a \code{listw} object created for example by \code{nb2listw}}

\item{zero.policy}{default \code{attr(listw, "zero.policy")} as set when \code{listw} was created, if attribute not set, use global option value; if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA}

\item{alternative}{a character string specifying the alternative hypothesis, must be one of \code{greater} (default), \code{less} or \code{two.sided}}

\item{na.action}{a function (default \code{na.fail}), can also be \code{na.omit} or \code{na.exclude} - in these cases the weights list will be subsetted to remove NAs in the data. It may be necessary to set \code{zero.policy} to TRUE because this subsetting may create no-neighbour observations. Note that only weights lists created without using the \code{glist} argument to \code{nb2listw} may be subsetted. If \code{na.pass} is used, zero is substituted for NA values in calculating the spatial lag}

\item{spChk}{should the data vector names be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use \code{get.spChkOption()}}

\item{adjust.n}{default TRUE, if FALSE the number of observations is not adjusted for no-neighbour observations, if TRUE, the number of observations is adjusted}
}
\value{
A list with class \code{htest} containing the following components: \describe{
\item{\code{statistic}}{the value of the standard deviate of normalized Moran's I.}
\item{\code{p.value}}{the p-value of the test.}
\item{\code{estimate}}{the value of the observed normalized Moran's I, its expectation and variance under normality.}
\item{\code{alternative}}{a character string describing the alternative hypothesis.}
\item{\code{method}}{a character string giving the assumption used for calculating the standard deviate.}
\item{\code{data.name}}{a character string giving the name(s) of the data.}
}
}
\description{
Normalized "\eqn{I_6}" Moran's test for spatial autocorrelation proposed by Tillé et al. (2026) using a spatial weights matrix in weights list form. The test assumes asymptotic normality, and results may be checked against those of \code{normalized.moran.mc} permutations.
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
coords.OLD <- cbind(COL.OLD$X, COL.OLD$Y)
normalized.moran.test(COL.OLD$CRIME, nb2listw(COL.nb, style="W"))
crime <- COL.OLD$CRIME
is.na(crime) <- sample(1:length(crime), 10)
res <- try(normalized.moran.test(crime, nb2listw(COL.nb, style="W"),
                                 na.action=na.fail))
normalized.moran.test(crime, nb2listw(COL.nb, style="W"), zero.policy=TRUE,
                      na.action=na.omit)
normalized.moran.test(crime, nb2listw(COL.nb, style="W"), zero.policy=TRUE,
                      na.action=na.exclude)
normalized.moran.test(crime, nb2listw(COL.nb, style="W"), na.action=na.pass)
}
\seealso{
\link{normalized.moran}, \link{normalized.moran.mc}
}
