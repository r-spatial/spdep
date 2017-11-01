\name{lm.morantest.sad}
\alias{lm.morantest.sad}
\alias{print.moransad}
\alias{summary.moransad}
\alias{print.summary.moransad}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Saddlepoint approximation of global Moran's I test}
\description{
  The function implements Tiefelsdorf's application of the Saddlepoint
approximation to global Moran's I's reference distribution. 
}
\usage{
lm.morantest.sad(model, listw, zero.policy=NULL, alternative="greater", 
  spChk=NULL, resfun=weighted.residuals, tol=.Machine$double.eps^0.5,
  maxiter=1000, tol.bounds=0.0001, zero.tol = 1e-07, Omega=NULL,
  save.M=NULL, save.U=NULL)
\method{print}{moransad}(x, ...)
\method{summary}{moransad}(object, ...)
\method{print}{summary.moransad}(x, ...) 
}
%- maybe also `usage' for other objects documented here.
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
  \item{tol}{the desired accuracy (convergence tolerance) for \code{uniroot}}
  \item{maxiter}{the maximum number of iterations for \code{uniroot}}
  \item{tol.bounds}{offset from bounds for \code{uniroot}}
  \item{zero.tol}{tolerance used to find eigenvalues close to absolute zero}
  \item{Omega}{A SAR process matrix may be passed in to test an alternative hypothesis, for example \code{Omega <- invIrW(listw, rho=0.1); Omega <- tcrossprod(Omega)}, \code{chol()} is taken internally}
  \item{save.M}{return the full M matrix for use in \code{spdep:::exactMoranAlt}}
  \item{save.U}{return the full U matrix for use in \code{spdep:::exactMoranAlt}}
  \item{x}{object to be printed}
  \item{object}{object to be summarised}
  \item{...}{arguments to be passed through}
}
\details{
  The function involves finding the eigenvalues of an n by n matrix, and
numerically finding the root for the Saddlepoint approximation, and should
therefore only be used with care when n is large.
}
\value{
  A list of class \code{moransad} with the following components:
  \item{statistic}{the value of the saddlepoint approximation of the
 standard deviate of global Moran's I.}
  \item{p.value}{the p-value of the test.}
  \item{estimate}{the value of the observed global Moran's I.}
  \item{alternative}{a character string describing the alternative hypothesis.}
  \item{method}{a character string giving the method used.}
  \item{data.name}{a character string giving the name(s) of the data.}
  \item{internal1}{Saddlepoint omega, r and u}
  \item{internal2}{f.root, iter and estim.prec from \code{uniroot}}
  \item{df}{degrees of freedom}
  \item{tau}{eigenvalues (excluding zero values)}
}

\references{Tiefelsdorf, M. 2002 The Saddlepoint approximation of Moran's I 
and local Moran's Ii reference distributions and their numerical evaluation. 
Geographical Analysis, 34, pp. 187--206.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{lm.morantest}}}

\examples{
if (require(rgdal, quietly=TRUE)) {
eire <- readOGR(system.file("shapes/eire.shp", package="spData")[1])
row.names(eire) <- as.character(eire$names)
proj4string(eire) <- CRS("+proj=utm +zone=30 +ellps=airy +units=km")
eire.nb <- poly2nb(eire)
e.lm <- lm(OWNCONS ~ ROADACC, data=eire)
lm.morantest(e.lm, nb2listw(eire.nb))
lm.morantest.sad(e.lm, nb2listw(eire.nb))
summary(lm.morantest.sad(e.lm, nb2listw(eire.nb)))
e.wlm <- lm(OWNCONS ~ ROADACC, data=eire, weights=RETSALE)
lm.morantest(e.wlm, nb2listw(eire.nb), resfun=rstudent)
lm.morantest.sad(e.wlm, nb2listw(eire.nb), resfun=rstudent)
}
}
\keyword{spatial}
