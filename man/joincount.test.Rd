% Copyright 2001-8 by Roger S. Bivand
\name{joincount.test}
\alias{joincount.test}
\alias{print.jclist}
\title{BB join count statistic for k-coloured factors}
\description{
The BB join count test for spatial autocorrelation using a spatial
weights matrix in weights list form for testing whether same-colour joins
occur more frequently than would be expected if the zones were labelled
in a spatially random way. The assumptions underlying the test are
sensitive to the form of the graph of neighbour relationships and other
factors, and results may be checked against those of \code{joincount.mc}
permutations.
}
\usage{
joincount.test(fx, listw, zero.policy=NULL, alternative="greater", spChk=NULL,
 adjust.n=TRUE)
\method{print}{jclist}(x, ...)
}
\arguments{
  \item{fx}{a factor of the same length as the neighbours and weights objects in listw}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{zero.policy}{default NULL, use global option value; if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA}
  \item{alternative}{a character string specifying the alternative hypothesis, must be one of "greater" (default), "less" or "two.sided".}
  \item{adjust.n}{default TRUE, if FALSE the number of observations is not adjusted for no-neighbour observations, if TRUE, the number of observations is adjusted consistently (up to and including spdep 0.3-28 the adjustment was inconsistent - thanks to Tomoki NAKAYA for a careful bug report)}
  \item{spChk}{should the data vector names be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use \code{get.spChkOption()}}
  \item{x}{object to be printed}
  \item{...}{arguments to be passed through for printing}
}

\value{
A list with class \code{jclist} of lists with class \code{htest} for each of the k colours containing the following components:
  \item{statistic}{the value of the standard deviate of the join count statistic.}
  \item{p.value}{the p-value of the test.}
  \item{estimate}{the value of the observed statistic, its expectation and variance under non-free sampling.}
  \item{alternative}{a character string describing the alternative hypothesis.}
  \item{method}{a character string giving the method used.}
  \item{data.name}{a character string giving the name(s) of the data.}
}
\references{Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 20.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}
\note{The derivation of the test (Cliff and Ord, 1981, p. 18) assumes that the weights matrix is symmetric. For inherently non-symmetric matrices, such as k-nearest neighbour matrices, \code{listw2U()} can be used to make the matrix symmetric. In non-symmetric weights matrix cases, the variance of the test statistic may be negative.
}

\seealso{\code{\link{joincount.mc}}, \code{\link{joincount.multi}}, \code{\link{listw2U}}}

\examples{
data(oldcol)
HICRIME <- cut(COL.OLD$CRIME, breaks=c(0,35,80), labels=c("low","high"))
names(HICRIME) <- rownames(COL.OLD)
joincount.test(HICRIME, nb2listw(COL.nb, style="B"))
joincount.test(HICRIME, nb2listw(COL.nb, style="C"))
joincount.test(HICRIME, nb2listw(COL.nb, style="S"))
joincount.test(HICRIME, nb2listw(COL.nb, style="W"))
by(card(COL.nb), HICRIME, summary)
print(is.symmetric.nb(COL.nb))
coords.OLD <- cbind(COL.OLD$X, COL.OLD$Y)
COL.k4.nb <- knn2nb(knearneigh(coords.OLD, 4))
print(is.symmetric.nb(COL.k4.nb))
joincount.test(HICRIME, nb2listw(COL.k4.nb, style="B"))
cat("Note non-symmetric weights matrix - use listw2U()\n")
joincount.test(HICRIME, listw2U(nb2listw(COL.k4.nb, style="B")))
}
\keyword{spatial}
