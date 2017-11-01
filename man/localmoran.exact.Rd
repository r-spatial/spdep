\name{localmoran.exact}
\alias{localmoran.exact}
\alias{localmoran.exact.alt}
\alias{print.localmoranex}
\alias{as.data.frame.localmoranex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Exact local Moran's Ii tests}
\description{
 \code{localmoran.exact} provides exact local Moran's Ii tests under the null hypothesis, while \code{localmoran.exact.alt} provides exact local Moran's Ii tests under the alternative hypothesis. In this case, the model may be a fitted model derived from a model fitted by \code{errorsarlm}, with the covariance matrix can be passed through the \code{Omega=} argument.
}
\usage{
localmoran.exact(model, select, nb, glist = NULL, style = "W", 
 zero.policy = NULL, alternative = "greater", spChk = NULL, 
 resfun = weighted.residuals, save.Vi = FALSE, useTP=FALSE, truncErr=1e-6, 
 zeroTreat=0.1)
localmoran.exact.alt(model, select, nb, glist = NULL, style = "W",
 zero.policy = NULL, alternative = "greater", spChk = NULL,
 resfun = weighted.residuals, Omega = NULL, save.Vi = FALSE,
 save.M = FALSE, useTP=FALSE, truncErr=1e-6, zeroTreat=0.1)
\method{print}{localmoranex}(x, \dots)
\method{as.data.frame}{localmoranex}(x, row.names=NULL, optional=FALSE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{an object of class \code{lm} returned by \code{lm} (assuming no global spatial autocorrelation), or an object of class \code{sarlm} returned by a spatial simultaneous autoregressive model fit (assuming global spatial autocorrelation represented by the model spatial coefficient); weights may be specified in the \code{lm} fit, but offsets should not be used}
  \item{select}{an integer vector of the id. numbers of zones to be tested; if missing, all zones}
  \item{nb}{a list of neighbours of class \code{nb}}
  \item{glist}{a list of general weights corresponding to neighbours}
  \item{style}{can take values W, B, C, and S}
  \item{zero.policy}{default NULL, use global option value; if TRUE assign zero to the lagged value of zones without
neighbours, if FALSE assign NA}
  \item{alternative}{a character string specifying the alternative hypothesis,
must be one of greater (default), less or two.sided.}
  \item{spChk}{should the data vector names be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use \code{get.spChkOption()}}
  \item{resfun}{default: weighted.residuals; the function to be used to extract residuals from the \code{lm} object, may be \code{residuals}, \code{weighted.residuals}, \code{rstandard}, or \code{rstudent}}
  \item{Omega}{A SAR process matrix may be passed in to test an alternative hypothesis, for example \code{Omega <- invIrW(listw, rho=0.1); Omega <- tcrossprod(Omega)}, \code{chol()} is taken internally}
  \item{save.Vi}{if TRUE, return the star-shaped weights lists for each zone 
tested}
  \item{save.M}{if TRUE, save a list of left and right M products}
  \item{useTP}{default FALSE, if TRUE, use truncation point in integration rather than upper=Inf, see Tiefelsdorf (2000), eq. 6.7, p.69}
  \item{truncErr}{when useTP=TRUE, pass truncation error to truncation point function}
  \item{zeroTreat}{when useTP=TRUE, pass zero adjustment to truncation point function}
  \item{x}{object to be printed}
  \item{row.names}{ignored argument to \code{as.data.frame.localmoranex};
row names assigned from localmoranex object}
  \item{optional}{ignored argument to \code{as.data.frame.localmoranex};
row names assigned from localmoranex object}
  \item{\dots}{arguments to be passed through}
}

\value{
A list with class \code{localmoranex} containing "select" lists, each with
 class \code{moranex} with the following components:
  \item{statistic}{the value of the exact standard deviate of global Moran's I.}
  \item{p.value}{the p-value of the test.}
  \item{estimate}{the value of the observed local Moran's Ii.}
  \item{method}{a character string giving the method used.}
  \item{alternative}{a character string describing the alternative hypothesis.}
  \item{gamma}{eigenvalues (two extreme values for null, vector for alternative)}
  \item{oType}{usually set to "E", but set to "N" if the integration leads to an out of domain value for \code{qnorm}, when the Normal assumption is substituted. This only occurs when the output p-value would be very close to zero}
  \item{data.name}{a character string giving the name(s) of the data.}
  \item{df}{degrees of freedom}
  \item{i}{zone tested}
  \item{Vi}{zone tested}

When the alternative is being tested, a list of left and right M products in attribute M.
}
\author{Markus Reder and Roger Bivand}

\seealso{\code{\link{lm.morantest.exact}}, \code{\link{localmoran.sad}}}
\examples{
if (require(rgdal, quietly=TRUE)) {
eire <- readOGR(system.file("shapes/eire.shp", package="spData")[1])
row.names(eire) <- as.character(eire$names)
proj4string(eire) <- CRS("+proj=utm +zone=30 +ellps=airy +units=km")
eire.nb <- poly2nb(eire)
e.lm <- lm(OWNCONS ~ ROADACC, data=eire)
localmoran.sad(e.lm, nb=eire.nb)
localmoran.exact(e.lm, nb=eire.nb)
localmoran.exact(e.lm, nb=eire.nb, useTP=TRUE)
e.errorsar <- errorsarlm(OWNCONS ~ ROADACC, data=eire,
 listw=nb2listw(eire.nb))
lm.target <- lm(e.errorsar$tary ~ e.errorsar$tarX - 1)
localmoran.exact.alt(lm.target, nb=eire.nb)
Omega <- invIrW(nb2listw(eire.nb), rho=0.6)
Omega1 <- tcrossprod(Omega)
localmoran.exact.alt(lm.target, nb=eire.nb, Omega=Omega1)
localmoran.exact.alt(lm.target, nb=eire.nb, Omega=Omega1, useTP=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}

