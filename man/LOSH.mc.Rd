\name{LOSH.mc}
\alias{LOSH.mc}
\title{Bootstrapping-based test for local spatial heteroscedasticity}
\description{
The function draws inferences about local spatial heteroscedasticity (LOSH) by means of the randomisation-based Monte-Carlo bootstrap proposed by Xu et al. (2014).
}
\usage{
LOSH.mc(x, listw, a = 2, nsim = 99, zero.policy = attr(listw, "zero.policy"),
 na.action = na.fail, spChk = NULL, adjust.n = TRUE, p.adjust.method = "none")
}
\arguments{
  \item{x}{a numeric vector of the same length as the neighbours list in listw}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{a}{the exponent applied to the local residuals; the default value of 2 leads to a measure of heterogeneity in the spatial variance}
  \item{nsim}{the number of randomisations used in the bootstrap}
  \item{zero.policy}{default \code{attr(listw, "zero.policy")} as set when \code{listw} was created, if attribute not set, use global option value; if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA}
  \item{na.action}{a function (default \code{na.fail}), can also be \code{na.omit} or \code{na.exclude} - in these cases the weights list will be subsetted to remove NAs in the data. It may be necessary to set zero.policy to TRUE because this subsetting may create no-neighbour observations. Note that only weights lists created without using the glist argument to \code{nb2listw} may be subsetted. If \code{na.pass} is used, zero is substituted for NA values in calculating the spatial lag. (Note that na.exclude will only work properly starting from R 1.9.0, na.omit and na.exclude assign the wrong classes in 1.8.*)}
  \item{spChk}{should the data vector names be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use \code{get.spChkOption()}}
  \item{adjust.n}{ default TRUE, if FALSE the number of observations is not adjusted for no-neighbour observations, if TRUE, the number of observations is adjusted}
  \item{p.adjust.method}{a character string specifying the probability value adjustment for multiple tests, default "none"; see \code{\link{p.adjustSP}}. Note that the number of multiple tests for each region is only taken as the number of neighbours + 1 for each region, rather than the total number of regions.}
}
\details{
The test calculates LOSH (see \code{\link{LOSH}}) and estimates pseudo p-values from a conditional bootstrap. Thereby, the i-th value in each location is held fixed, whereas all other values are permuted \code{nsim} times over all other spatial units.
}
  \value{
    \item{Hi}{LOSH statistic}
    \item{E.Hi}{expectation of LOSH}
    \item{Var.Hi}{variance of LOSH}
    \item{Z.Hi}{the approximately chi-square distributed test statistics}
    \item{x_bar_i}{local spatially weighted mean values}
    \item{ei}{residuals about local spatially weighted mean values}
    \item{Pr()}{p-values for \code{Hi} obtained from a conditional bootstrap distribution}
  }
\references{
Ord, J. K., & Getis, A. 2012. Local spatial heteroscedasticity (LOSH), The Annals of Regional Science, 48 (2), 529--539; Xu, M., Mei, C. L., & Yan, N. 2014. A note on the null distribution of the local spatial heteroscedasticity (LOSH) statistic. The Annals of Regional Science, 52 (3), 697--710.
}

\author{René Westerholt \email{rene.westerholt@tu-dortmund.de}}

\seealso{\code{\link{LOSH}}, \code{\link{LOSH.mc}}}
  
\examples{
    data(columbus, package="spData")
    resLOSH_mc <- LOSH.mc(columbus$CRIME, nb2listw(col.gal.nb), 2, 100)
    summary(resLOSH_mc)
    resLOSH_cs <- LOSH.cs(columbus$CRIME, nb2listw(col.gal.nb))
    summary(resLOSH_cs)
    plot(resLOSH_mc[,"Pr()"], resLOSH_cs[,"Pr()"])
}
\keyword{spatial}
