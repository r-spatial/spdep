% Copyright 2002-4 by Roger S. Bivand
\name{localmoran.sad}
\alias{localmoran.sad}
\alias{listw2star}
\alias{print.summary.localmoransad}
\alias{summary.localmoransad}
\alias{print.localmoransad}
\alias{as.data.frame.localmoransad}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Saddlepoint approximation of local Moran's Ii tests}
\description{
  The function implements Tiefelsdorf's application of the Saddlepoint
approximation to local Moran's Ii's reference distribution. If the model
object is of class "lm", global independence is assumed; if of class
"sarlm", global dependence is assumed to be represented by the spatial
parameter of that model. Tests are reported separately for each zone
selected, and may be summarised using \code{summary.localmoransad}.
Values of local Moran's Ii agree with those from \code{localmoran()},
but in that function, the standard deviate - here the Saddlepoint
approximation - is based on the randomisation assumption.
}
\usage{
localmoran.sad(model, select, nb, glist=NULL, style="W",
 zero.policy=NULL, alternative="greater", spChk=NULL,
 resfun=weighted.residuals, save.Vi=FALSE,
 tol = .Machine$double.eps^0.5, maxiter = 1000, tol.bounds=0.0001,
 save.M=FALSE, Omega = NULL)
%as.data.frame.localmoransad(x, row.names=NULL, optional=FALSE) 
\method{print}{localmoransad}(x, ...)
\method{summary}{localmoransad}(object, ...)
\method{print}{summary.localmoransad}(x, ...)
listw2star(listw, ireg, style, n, D, a, zero.policy=NULL)
}
%- maybe also `usage' for other objects documented here.
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
  \item{save.Vi}{if TRUE, return the star-shaped weights lists for each zone 
tested}
  \item{tol}{the desired accuracy (convergence tolerance) for \code{uniroot}}
  \item{maxiter}{the maximum number of iterations for \code{uniroot}}
  \item{tol.bounds}{offset from bounds for \code{uniroot}}
  \item{save.M}{if TRUE, save a list of left and right M products in a list for the conditional tests, or a list of the regression model matrix components}
  \item{Omega}{A SAR process matrix may be passed in to test an alternative hypothesis, for example \code{Omega <- invIrW(listw, rho=0.1); Omega <- tcrossprod(Omega)}, \code{chol()} is taken internally}
  \item{x}{object to be printed}
%  \item{row.names}{ignored argument to \code{as.data.frame.localmoransad};
%row names assigned from localmoransad object}
%  \item{optional}{ignored argument to \code{as.data.frame.localmoransad};
%row names assigned from localmoransad object}
  \item{object}{object to be summarised}
  \item{...}{arguments to be passed through}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{ireg}{a zone number}
  \item{n}{internal value depending on listw and style}
  \item{D}{internal value depending on listw and style}
  \item{a}{internal value depending on listw and style}
}
\details{
  The function implements the analytical eigenvalue calculation together with
trace shortcuts given or suggested in Tiefelsdorf (2002), partly
following remarks by J. Keith Ord, and uses the
Saddlepoint analytical solution from Tiefelsdorf's SPSS code. 

If a histogram of the probability values of the saddlepoint estimate
for the assumption of global independence is not approximately flat,
the assumption is probably unjustified, and re-estimation with global
dependence is recommended.

No n by n matrices are needed at any point for the test assuming no global
dependence, the star-shaped weights matrices being handled as listw
lists. When the test is made on residuals from a spatial regression,
taking a global process into account. n by n matrices are necessary,
and memory constraints may be reached for large lattices.
}
\value{
A list with class \code{localmoransad} containing "select" lists, each with
 class \code{moransad} with the following components:
  \item{statistic}{the value of the saddlepoint approximation of the
 standard deviate of local Moran's Ii.}
  \item{p.value}{the p-value of the test.}
  \item{estimate}{the value of the observed local Moran's Ii.}
  \item{alternative}{a character string describing the alternative hypothesis.}
  \item{method}{a character string giving the method used.}
  \item{data.name}{a character string giving the name(s) of the data.}
  \item{internal1}{Saddlepoint omega, r and u}
  \item{df}{degrees of freedom}
  \item{tau}{maximum and minimum analytical eigenvalues}
  \item{i}{zone tested}
}
\references{Tiefelsdorf, M. 2002 The Saddlepoint approximation of Moran's I 
and local Moran's Ii reference distributions and their numerical evaluation. 
Geographical Analysis, 34, pp. 187--206.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{localmoran}}, \code{\link{lm.morantest}},
 \code{\link{lm.morantest.sad}}, \code{\link{errorsarlm}}}

\examples{
if (require(rgdal, quietly=TRUE)) {
eire <- readOGR(system.file("shapes/eire.shp", package="spData")[1])
row.names(eire) <- as.character(eire$names)
proj4string(eire) <- CRS("+proj=utm +zone=30 +ellps=airy +units=km")
eire.nb <- poly2nb(eire)
lw <- nb2listw(eire.nb)
e.lm <- lm(OWNCONS ~ ROADACC, data=eire)
e.locmor <- summary(localmoran.sad(e.lm, nb=eire.nb))
e.locmor
mean(e.locmor[,1])
sum(e.locmor[,1])/Szero(lw)
lm.morantest(e.lm, lw)
# note equality for mean() only when the sum of weights equals
# the number of observations (thanks to Juergen Symanzik)
hist(e.locmor[,"Pr. (Sad)"])
e.wlm <- lm(OWNCONS ~ ROADACC, data=eire, weights=RETSALE)
e.locmorw1 <- summary(localmoran.sad(e.wlm, nb=eire.nb, resfun=weighted.residuals))
e.locmorw1
e.locmorw2 <- summary(localmoran.sad(e.wlm, nb=eire.nb, resfun=rstudent))
e.locmorw2
e.errorsar <- errorsarlm(OWNCONS ~ ROADACC, data=eire,
  listw=lw)
e.errorsar
lm.target <- lm(e.errorsar$tary ~ e.errorsar$tarX - 1)
Omega <- tcrossprod(invIrW(lw, rho=e.errorsar$lambda))
e.clocmor <- summary(localmoran.sad(lm.target, nb=eire.nb, Omega=Omega))
e.clocmor
hist(e.clocmor[,"Pr. (Sad)"])
}
}
\keyword{spatial}
