\name{LOSH.cs}
\alias{LOSH.cs}
\title{Chi-square based test for local spatial heteroscedasticity}
\description{
The function implements the chi-square based test statistic for local spatial heteroscedasticity (LOSH) as proposed by Ord & Getis (2012).
}
\usage{
LOSH.cs(x, listw, zero.policy = NULL, na.action = na.fail, 
                 p.adjust.method = "none", spChk = NULL)
}
\arguments{
  \item{x}{a numeric vector of the same length as the neighbours list in listw}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{zero.policy}{default NULL, use global option value; if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA}
  \item{na.action}{a function (default \code{na.fail}), can also be \code{na.omit} or \code{na.exclude} - in these cases the weights list will be subsetted to remove NAs in the data. It may be necessary to set zero.policy to TRUE because this subsetting may create no-neighbour observations. Note that only weights lists created without using the glist argument to \code{nb2listw} may be subsetted. If \code{na.pass} is used, zero is substituted for NA values in calculating the spatial lag. (Note that na.exclude will only work properly starting from R 1.9.0, na.omit and na.exclude assign the wrong classes in 1.8.*)}
 \item{p.adjust.method}{a character string specifying the probability value adjustment for multiple tests, default "none"; see \code{\link{p.adjustSP}}. Note that the number of multiple tests for each region is only taken as the number of neighbours + 1 for each region, rather than the total number of regions.}
  \item{spChk}{should the data vector names be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use \code{get.spChkOption()}}
}
\details{
The test uses \emph{a = 2} (see \code{\link{LOSH}}) because chi-square based inference is not applicable with other exponents. The function makes use of \code{\link{LOSH}} in its calculations.
}
  \value{
    \item{Hi}{LOSH statistic}
    \item{E.Hi}{expectation of LOSH}
    \item{Var.Hi}{variance of LOSH}
    \item{Z.Hi}{the approximately chi-square distributed test statistics}
    \item{x_bar_i}{local spatially weighted mean values}
    \item{ei}{residuals about local spatially weighted mean values}
    \item{Pr()}{p-values for \code{Hi} obtained from a non-central Chi-square distribution with \eqn{2/Var.Hi} degrees of freedom}
}
\references{
Ord, J. K., & Getis, A. 2012. Local spatial heteroscedasticity (LOSH), The Annals of Regional Science, 48 (2), 529--539.
}

\author{Rene Westerholt \email{westerholt@uni-heidelberg.de}}

\seealso{\code{\link{LOSH}}, \code{\link{LOSH.mc}}}
  
\examples{
    data(boston, package="spData")
    resLOSH <- LOSH.cs(boston.c$NOX, nb2listw(boston.soi))
    hist(resLOSH[,"Hi"])
    mean(resLOSH[,"Hi"])
}
\keyword{spatial}
