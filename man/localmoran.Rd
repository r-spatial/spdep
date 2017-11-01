% Copyright 2001-8 by Roger S. Bivand
\name{localmoran}
\alias{localmoran}
\title{Local Moran's I statistic}
\description{
 The local spatial statistic Moran's I is calculated for each zone based on the
spatial weights object used. The values returned include a Z-value, and may be
used as a diagnostic tool. The statistic is:
\deqn{I_i = \frac{(x_i-\bar{x})}{{\sum_{k=1}^{n}(x_k-\bar{x})^2}/(n-1)}{\sum_{j=1}^{n}w_{ij}(x_j-\bar{x})}},
and its expectation and variance are given in Anselin (1995).
}
\usage{
localmoran(x, listw, zero.policy=NULL, na.action=na.fail, 
	alternative = "greater", p.adjust.method="none", mlvar=TRUE,
        spChk=NULL, sokal98=FALSE)
}
\arguments{
  \item{x}{a numeric vector the same length as the neighbours list in listw}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{zero.policy}{default NULL, use global option value; if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA}
  \item{na.action}{a function (default \code{na.fail}), can also be \code{na.omit} or \code{na.exclude} - in these cases the weights list will be subsetted to remove NAs in the data. It may be necessary to set zero.policy to TRUE because this subsetting may create no-neighbour observations. Note that only weights lists created without using the glist argument to \code{nb2listw} may be subsetted. If \code{na.pass} is used, zero is substituted for NA values in calculating the spatial lag. (Note that na.exclude will only work properly starting from R 1.9.0, na.omit and na.exclude assign the wrong classes in 1.8.*)}
  \item{alternative}{a character string specifying the alternative hypothesis, must be one of greater (default), less or two.sided.}
  \item{p.adjust.method}{a character string specifying the probability value adjustment for multiple tests, default "none"; see \code{\link{p.adjustSP}}. Note that the number of multiple tests for each region is only taken as the number of neighbours + 1 for each region, rather than the total number of regions.}
  \item{mlvar}{default TRUE: values of local Moran's I are reported using the variance of the variable of interest (sum of squared deviances over n), but can be reported as the sample variance, dividing by (n-1) instead; both are used in other implementations.}
  \item{spChk}{should the data vector names be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use \code{get.spChkOption()}}
  \item{sokal98}{default=FALSE, if TRUE, use VI as defined by Sokal et al. (1998)}
}

\details{
The values of local Moran's I are divided by the variance (or sample variance) of the variable of interest to accord with Table 1, p. 103, and formula (12), p. 99, in Anselin (1995), rathar than his formula (7), p. 98.
}

\value{

  \item{Ii}{local moran statistic}
  \item{E.Ii}{expectation of local moran statistic}
  \item{Var.Ii}{variance of local moran statistic}
  \item{Z.Ii}{standard deviate of local moran statistic}
  \item{Pr()}{p-value of local moran statistic}
}
\references{Anselin, L. 1995. Local indicators of spatial association,
Geographical Analysis, 27, 93--115;
Getis, A. and Ord, J. K. 1996 Local spatial
statistics: an overview. In P. Longley and M. Batty (eds) \emph{Spatial
analysis: modelling in a GIS environment} (Cambridge: Geoinformation
International), 261--277; Sokal, R. R, Oden, N. L. and Thomson, B. A. 1998. Local Spatial Autocorrelation in a Biological Model. Geographical Analysis, 30. 331--354.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{localG}}}

\examples{
data(afcon, package="spData")
oid <- order(afcon$id)
resI <- localmoran(afcon$totcon, nb2listw(paper.nb))
printCoefmat(data.frame(resI[oid,], row.names=afcon$name[oid]),
 check.names=FALSE)
hist(resI[,5])
mean(resI[,1])
sum(resI[,1])/Szero(nb2listw(paper.nb))
moran.test(afcon$totcon, nb2listw(paper.nb))
# note equality for mean() only when the sum of weights equals
# the number of observations (thanks to Juergen Symanzik)
resI <- localmoran(afcon$totcon, nb2listw(paper.nb),
 p.adjust.method="bonferroni")
printCoefmat(data.frame(resI[oid,], row.names=afcon$name[oid]),
 check.names=FALSE)
hist(resI[,5])
totcon <-afcon$totcon
is.na(totcon) <- sample(1:length(totcon), 5)
totcon
resI.na <- localmoran(totcon, nb2listw(paper.nb), na.action=na.exclude,
 zero.policy=TRUE)
if (class(attr(resI.na, "na.action")) == "exclude") {
 print(data.frame(resI.na[oid,], row.names=afcon$name[oid]), digits=2)
} else print(resI.na, digits=2)
resG <- localG(afcon$totcon, nb2listw(include.self(paper.nb)))
print(data.frame(resG[oid], row.names=afcon$name[oid]), digits=2)

}
\keyword{spatial}
