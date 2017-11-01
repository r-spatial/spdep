\name{lm.morantest.exact}
\alias{lm.morantest.exact}
\alias{print.moranex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Exact global Moran's I test}
\description{
  The function implements Tiefelsdorf's exact global Moran's I test.
}
\usage{
lm.morantest.exact(model, listw, zero.policy = NULL, alternative = "greater",
 spChk = NULL, resfun = weighted.residuals, zero.tol = 1e-07, Omega=NULL,
 save.M=NULL, save.U=NULL, useTP=FALSE, truncErr=1e-6, zeroTreat=0.1)
\method{print}{moranex}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{an object of class \code{lm} returned by \code{lm}; weights
may be specified in the \code{lm} fit, but offsets should not be used}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{zero.policy}{default NULL, use global option value; if TRUE assign zero to the lagged value of zones without
neighbours, if FALSE assign NA}
  \item{alternative}{a character string specifying the alternative hypothesis,
must be one of greater (default), less or two.sided.}
  \item{spChk}{should the data vector names be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use \code{get.spChkOption()}}
  \item{resfun}{default: weighted.residuals; the function to be used to extract residuals from the \code{lm} object, may be \code{residuals}, \code{weighted.residuals}, \code{rstandard}, or \code{rstudent}}
  \item{zero.tol}{tolerance used to find eigenvalues close to absolute zero}
  \item{Omega}{A SAR process matrix may be passed in to test an alternative hypothesis, for example \code{Omega <- invIrW(listw, rho=0.1); Omega <- tcrossprod(Omega)}, \code{chol()} is taken internally}
  \item{save.M}{return the full M matrix for use in \code{spdep:::exactMoranAlt}}
  \item{save.U}{return the full U matrix for use in \code{spdep:::exactMoranAlt}}
  \item{useTP}{default FALSE, if TRUE, use truncation point in integration rather than upper=Inf, see Tiefelsdorf (2000), eq. 6.7, p.69}
  \item{truncErr}{when useTP=TRUE, pass truncation error to truncation point function}
  \item{zeroTreat}{when useTP=TRUE, pass zero adjustment to truncation point function}
  \item{x}{a moranex object}
  \item{\dots}{arguments to be passed through}
}
\value{
  A list of class \code{moranex} with the following components:
  \item{statistic}{the value of the saddlepoint approximation of the
 standard deviate of global Moran's I.}
  \item{p.value}{the p-value of the test.}
  \item{estimate}{the value of the observed global Moran's I.}
  \item{method}{a character string giving the method used.}
  \item{alternative}{a character string describing the alternative hypothesis.}
  \item{gamma}{eigenvalues (excluding zero values)}
  \item{oType}{usually set to "E"}
  \item{data.name}{a character string giving the name(s) of the data.}
  \item{df}{degrees of freedom}
}
\author{Markus Reder and Roger Bivand}

\references{Roger Bivand, Werner G. Müller and Markus Reder (2009) "Power calculations for global and local Moran's I." \emph{Computational Statistics \& Data Analysis} 53, 2859-2872.}

\seealso{\code{\link{lm.morantest.sad}}}
\examples{
if (require(rgdal, quietly=TRUE)) {
eire <- readOGR(system.file("shapes/eire.shp", package="spData")[1])
row.names(eire) <- as.character(eire$names)
proj4string(eire) <- CRS("+proj=utm +zone=30 +ellps=airy +units=km")
eire.nb <- poly2nb(eire)
e.lm <- lm(OWNCONS ~ ROADACC, data=eire)
lm.morantest(e.lm, nb2listw(eire.nb))
lm.morantest.sad(e.lm, nb2listw(eire.nb))
lm.morantest.exact(e.lm, nb2listw(eire.nb))
lm.morantest.exact(e.lm, nb2listw(eire.nb), useTP=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}

