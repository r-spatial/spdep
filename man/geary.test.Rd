% Copyright 2001-9 by Roger S. Bivand
\name{geary.test}
\alias{geary.test}
\title{Geary's C test for spatial autocorrelation}
\description{
   Geary's test for spatial autocorrelation using a spatial weights matrix in weights list form. The assumptions underlying the test are sensitive to the form of the graph of neighbour relationships and other factors, and results may be checked against those of \code{geary.mc} permutations. 
}
\usage{
geary.test(x, listw, randomisation=TRUE, zero.policy=NULL,
    alternative="greater", spChk=NULL, adjust.n=TRUE)
}

\arguments{
  \item{x}{a numeric vector the same length as the neighbours list in listw}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{randomisation}{variance of I calculated under the assumption of randomisation, if FALSE normality}
  \item{zero.policy}{default NULL, use global option value; if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA}
  \item{alternative}{a character string specifying the alternative hypothesis, must be one of "greater" (default), "less" or "two.sided".}
  \item{spChk}{should the data vector names be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use \code{get.spChkOption()}}
  \item{adjust.n}{default TRUE, if FALSE the number of observations is not adjusted for no-neighbour observations, if TRUE, the number of observations is adjusted}
}

\value{
A list with class \code{htest} containing the following components:
  \item{statistic}{the value of the standard deviate of Geary's C, in the order given in Cliff and Ord 1973, p. 21, which is (EC - C) / sqrt(VC), that is with the sign reversed with respect to the more usual (C - EC) / sqrt(VC); this means that the \dQuote{greater} alternative for the Geary C test corresponds to the \dQuote{greater} alternative for Moran's I test.}
  \item{p.value}{the p-value of the test.}
  \item{estimate}{the value of the observed Geary's C, its expectation and 
 variance under the method assumption.}
  \item{alternative}{a character string describing the alternative hypothesis.}
  \item{method}{a character string giving the assumption used for calculating 
 the standard deviate.}
  \item{data.name}{a character string giving the name(s) of the data.}
}
\references{Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 21, Cliff, A. D., Ord, J. K. 1973 Spatial Autocorrelation, Pion, pp. 15-16, 21.}

\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\note{The derivation of the test (Cliff and Ord, 1981, p. 18) assumes that the weights matrix is symmetric. For inherently non-symmetric matrices, such as k-nearest neighbour matrices, \code{listw2U()} can be used to make the matrix symmetric. In non-symmetric weights matrix cases, the variance of the test statistic may be negative (thanks to Franz Munoz I for a well documented bug report). Geary's C is affected by non-symmetric weights under normality much more than Moran's I. From 0.4-35, the sign of the standard deviate of C is changed to match Cliff and Ord (1973, p. 21).
}

\seealso{\code{\link{geary}}, \code{\link{geary.mc}}, \code{\link{listw2U}}}

\examples{
data(oldcol)
geary.test(COL.OLD$CRIME, nb2listw(COL.nb, style="W"))
geary.test(COL.OLD$CRIME, nb2listw(COL.nb, style="W"),
 randomisation=FALSE)
colold.lags <- nblag(COL.nb, 3)
geary.test(COL.OLD$CRIME, nb2listw(colold.lags[[2]],
 style="W"))
geary.test(COL.OLD$CRIME, nb2listw(colold.lags[[3]],
 style="W"), alternative="greater")
print(is.symmetric.nb(COL.nb))
coords.OLD <- cbind(COL.OLD$X, COL.OLD$Y)
COL.k4.nb <- knn2nb(knearneigh(coords.OLD, 4))
print(is.symmetric.nb(COL.k4.nb))
geary.test(COL.OLD$CRIME, nb2listw(COL.k4.nb, style="W"))
geary.test(COL.OLD$CRIME, nb2listw(COL.k4.nb, style="W"),
 randomisation=FALSE)
cat("Note non-symmetric weights matrix - use listw2U()\n")
geary.test(COL.OLD$CRIME, listw2U(nb2listw(COL.k4.nb,
 style="W")))
geary.test(COL.OLD$CRIME, listw2U(nb2listw(COL.k4.nb,
 style="W")), randomisation=FALSE)
}
\keyword{spatial}
