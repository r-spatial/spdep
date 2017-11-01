% Copyright 2001-8 by Roger S. Bivand
\name{moran.test}
\alias{moran.test}
\title{Moran's I test for spatial autocorrelation}
\description{
 Moran's test for spatial autocorrelation using a spatial weights matrix in weights list form. The assumptions underlying the test are sensitive to the form of the graph of neighbour relationships and other factors, and results may be checked against those of \code{moran.mc} permutations. 
}
\usage{
moran.test(x, listw, randomisation=TRUE, zero.policy=NULL,
 alternative="greater", rank = FALSE, na.action=na.fail, spChk=NULL, adjust.n=TRUE)
}
\arguments{
  \item{x}{a numeric vector the same length as the neighbours list in listw}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{randomisation}{variance of I calculated under the assumption of randomisation, if FALSE normality}
  \item{zero.policy}{default NULL, use global option value; if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA}
  \item{alternative}{a character string specifying the alternative hypothesis, must be one of greater (default), less or two.sided.}
  \item{rank}{logical value - default FALSE for continuous variables, if TRUE, uses the adaptation of Moran's I for ranks suggested by Cliff and Ord (1981, p. 46)}
  \item{na.action}{a function (default \code{na.fail}), can also be \code{na.omit} or \code{na.exclude} - in these cases the weights list will be subsetted to remove NAs in the data. It may be necessary to set zero.policy to TRUE because this subsetting may create no-neighbour observations. Note that only weights lists created without using the glist argument to \code{nb2listw} may be subsetted. If \code{na.pass} is used, zero is substituted for NA values in calculating the spatial lag}
  \item{spChk}{should the data vector names be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use \code{get.spChkOption()}}
  \item{adjust.n}{default TRUE, if FALSE the number of observations is not adjusted for no-neighbour observations, if TRUE, the number of observations is adjusted}
}


\value{
A list with class \code{htest} containing the following components:
  \item{statistic}{the value of the standard deviate of Moran's I.}
  \item{p.value}{the p-value of the test.}
  \item{estimate}{the value of the observed Moran's I, its expectation and variance under the method assumption.}
  \item{alternative}{a character string describing the alternative hypothesis.}
  \item{method}{a character string giving the assumption used for calculating the standard deviate.}
  \item{data.name}{a character string giving the name(s) of the data.}
}
\note{Var(I) is taken from Cliff and Ord (1969, p. 28), and 
Goodchild's CATMOG 47 (1986),
see also Upton & Fingleton (1985) p. 171; it agrees with SpaceStat,
see Tutorial workbook Chapter 22; VI is the second crude moment minus the
square of the first crude moment. The derivation of the test (Cliff and Ord, 1981, p. 18) assumes that the weights matrix is symmetric. For inherently non-symmetric matrices, such as k-nearest neighbour matrices, \code{listw2U()} can be used to make the matrix symmetric.
}

\references{Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 21.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{moran}}, \code{\link{moran.mc}}, \code{\link{listw2U}}}

\examples{
data(oldcol)
coords.OLD <- cbind(COL.OLD$X, COL.OLD$Y)
moran.test(COL.OLD$CRIME, nb2listw(COL.nb, style="W"))
moran.test(COL.OLD$CRIME, nb2listw(COL.nb, style="B"))
moran.test(COL.OLD$CRIME, nb2listw(COL.nb, style="C"))
moran.test(COL.OLD$CRIME, nb2listw(COL.nb, style="S"))
moran.test(COL.OLD$CRIME, nb2listw(COL.nb, style="W"),
 randomisation=FALSE)
colold.lags <- nblag(COL.nb, 3)
moran.test(COL.OLD$CRIME, nb2listw(colold.lags[[2]],
 style="W"))
moran.test(COL.OLD$CRIME, nb2listw(colold.lags[[3]],
 style="W"))
print(is.symmetric.nb(COL.nb))
COL.k4.nb <- knn2nb(knearneigh(coords.OLD, 4))
print(is.symmetric.nb(COL.k4.nb))
moran.test(COL.OLD$CRIME, nb2listw(COL.k4.nb, style="W"))
moran.test(COL.OLD$CRIME, nb2listw(COL.k4.nb, style="W"),
 randomisation=FALSE)
cat("Note: non-symmetric weights matrix, use listw2U()")
moran.test(COL.OLD$CRIME, listw2U(nb2listw(COL.k4.nb,
 style="W")))
moran.test(COL.OLD$CRIME, listw2U(nb2listw(COL.k4.nb,
 style="W")), randomisation=FALSE)
ranks <- rank(COL.OLD$CRIME)
names(ranks) <- rownames(COL.OLD)
moran.test(ranks, nb2listw(COL.nb, style="W"), rank=TRUE)
crime <- COL.OLD$CRIME
is.na(crime) <- sample(1:length(crime), 10)
res <- try(moran.test(crime, nb2listw(COL.nb, style="W"),
 na.action=na.fail))
res
moran.test(crime, nb2listw(COL.nb, style="W"), zero.policy=TRUE,
 na.action=na.omit)
moran.test(crime, nb2listw(COL.nb, style="W"), zero.policy=TRUE,
 na.action=na.exclude)
moran.test(crime, nb2listw(COL.nb, style="W"), na.action=na.pass)
}
\keyword{spatial}
