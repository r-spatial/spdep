% Copyright 2001-8 by Roger S. Bivand
\name{localmoran}
\alias{localmoran}
\alias{localmoran_perm}
\title{Local Moran's I statistic}
\description{
 The local spatial statistic Moran's I is calculated for each zone based on the
spatial weights object used. The values returned include a Z-value, and may be
used as a diagnostic tool. The statistic and its expectation and variance were given in Anselin (1995), but those from Sokal et al. (1998) are implemented here.
}
\usage{
localmoran(x, listw, zero.policy=attr(listw, "zero.policy"), na.action=na.fail,
        conditional=TRUE, alternative = "two.sided", mlvar=TRUE,
        spChk=NULL, adjust.x=FALSE)
localmoran_perm(x, listw, nsim=499, zero.policy=attr(listw, "zero.policy"), 
        na.action=na.fail, alternative = "two.sided", mlvar=TRUE,
        spChk=NULL, adjust.x=FALSE, sample_Ei=TRUE, iseed=NULL,
        no_repeat_in_row=FALSE)
}
\arguments{
  \item{x}{a numeric vector the same length as the neighbours list in listw}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{zero.policy}{default default \code{attr(listw, "zero.policy")} as set when \code{listw} was created, if attribute not set, use global option value; if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA}
  \item{na.action}{a function (default \code{na.fail}), can also be \code{na.omit} or \code{na.exclude} - in these cases the weights list will be subsetted to remove NAs in the data. It may be necessary to set zero.policy to TRUE because this subsetting may create no-neighbour observations. Note that only weights lists created without using the glist argument to \code{nb2listw} may be subsetted. If \code{na.pass} is used, zero is substituted for NA values in calculating the spatial lag. (Note that na.exclude will only work properly starting from R 1.9.0, na.omit and na.exclude assign the wrong classes in 1.8.*)}
  \item{conditional}{default TRUE: expectation and variance are calculated using the conditional randomization null (Sokal 1998 Eqs. A7 & A8). Elaboration of these changes available in Sauer et al. (2021). If FALSE: expectation and variance are calculated using the total randomization null (Sokal 1998 Eqs. A3 & A4).}
  \item{alternative}{a character string specifying the alternative hypothesis, must be one of greater, less or two.sided (default).}
  \item{mlvar}{default TRUE: values of local Moran's I are reported using the variance of the variable of interest (sum of squared deviances over n), but can be reported as the sample variance, dividing by (n-1) instead; both are used in other implementations.}
  \item{spChk}{should the data vector names be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use \code{get.spChkOption()}}
  \item{adjust.x}{default FALSE, if TRUE, x values of observations with no neighbours are omitted in the mean of x}
  \item{nsim}{default 499, number of conditonal permutation simulations}
  \item{sample_Ei}{default TRUE; if conditional permutation, use the sample $E_i$ values, or the analytical values, leaving only variances calculated by simulation.}
  \item{iseed}{default NULL, used to set the seed; the output will only be reproducible if the count of CPU cores across which computation is distributed is the same}
  \item{no_repeat_in_row}{default \code{FALSE}, if \code{TRUE}, sample conditionally in each row without replacements to avoid duplicate values, \url{https://github.com/r-spatial/spdep/issues/124}}
}

\details{
The values of local Moran's I are divided by the variance (or sample variance) of the variable of interest to accord with Table 1, p. 103, and formula (12), p. 99, in Anselin (1995), rather than his formula (7), p. 98. The variance of the local Moran statistic is taken from Sokal et al. (1998) p. 334, equations 4 & 5 or equations 7 & 8 located depending on user specification. By default, the implementation divides by n, not (n-1) in calculating the variance and higher moments. Conditional code contributed by Jeff Sauer and Levi Wolf.
}

\value{

  \item{Ii}{local moran statistic}
  \item{E.Ii}{expectation of local moran statistic; for \code{localmoran_perm}the permutation sample means}
  \item{Var.Ii}{variance of local moran statistic; for \code{localmoran_perm}the permutation sample standard deviations}
  \item{Z.Ii}{standard deviate of local moran statistic; for \code{localmoran_perm} based on permutation sample means and standard deviations}
  \item{Pr()}{p-value of local moran statistic using \code{pnorm()}; for \code{localmoran_perm} using standard deviatse based on permutation sample means and standard deviations}
  \item{Pr() Sim}{For \code{localmoran_perm}, \code{rank()} and \code{punif()} of observed statistic rank for [0, 1] p-values using \code{alternative=}}
  \item{Pr(folded) Sim}{the simulation folded [0, 0.5] range ranked p-value (based on \url{https://github.com/pysal/esda/blob/4a63e0b5df1e754b17b5f1205b8cadcbecc5e061/esda/crand.py#L211-L213})}
  \item{Skewness}{For \code{localmoran_perm}, the output of \code{e1071::skewness()} for the permutation samples underlying the standard deviates}
  \item{Kurtosis}{For \code{localmoran_perm}, the output of \code{e1071::kurtosis()} for the permutation samples underlying the standard deviates}
In addition, an attribute data frame \code{"quadr"} with mean and median quadrant columns, and a column splitting on the demeaned variable and lagged demeaned variable at zero.
}
\note{Conditional permutations added for comparative purposes; permutations are over the whole data vector omitting the observation itself. For p-value adjustment, use \code{p.adjust()} or \code{p.adjustSP()} on the output vector.}
\references{Anselin, L. 1995. Local indicators of spatial association,
Geographical Analysis, 27, 93--115;
Getis, A. and Ord, J. K. 1996 Local spatial
statistics: an overview. In P. Longley and M. Batty (eds) \emph{Spatial
analysis: modelling in a GIS environment} (Cambridge: Geoinformation
International), 261--277; Sokal, R. R, Oden, N. L. and Thomson, B. A. 1998. Local Spatial Autocorrelation in a Biological Model. Geographical Analysis, 30. 331--354; 
Bivand RS, Wong DWS 2018 Comparing implementations of global and local indicators of spatial association. TEST, 27(3), 716--748 \doi{10.1007/s11749-018-0599-x}; 
Sauer, J., Oshan, T. M., Rey, S., & Wolf, L. J. 2021. The Importance of Null Hypotheses: Understanding Differences in Local Moran’s under Heteroskedasticity. Geographical Analysis. \doi{10.1111/gean.12304}

Bivand, R. (2022), R Packages for Analyzing Spatial Data: A Comparative Case Study with Areal Data. Geographical Analysis, 54(3), 488-518. \doi{10.1111/gean.12319}

}
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
resI <- localmoran(afcon$totcon, nb2listw(paper.nb))
printCoefmat(data.frame(resI[oid,], row.names=afcon$name[oid]),
 check.names=FALSE)
hist(p.adjust(resI[,5], method="bonferroni"))
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
set.seed(1)
resI_p <- localmoran_perm(afcon$totcon, nb2listw(paper.nb))
printCoefmat(data.frame(resI_p[oid,], row.names=afcon$name[oid]),
 check.names=FALSE)
}
\keyword{spatial}
