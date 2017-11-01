\name{GMerrorsar}
\alias{GMerrorsar}
\alias{residuals.gmsar}
\alias{deviance.gmsar}
\alias{coef.gmsar}
\alias{fitted.gmsar}
\alias{print.gmsar}
\alias{summary.gmsar}
\alias{print.summary.gmsar}
\alias{GMargminImage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Spatial simultaneous autoregressive error model estimation by GMM}
\description{
  An implementation of Kelejian and Prucha's generalised moments estimator for the autoregressive parameter in a spatial model.
}
\usage{
GMerrorsar(formula, data = list(), listw, na.action = na.fail,
 zero.policy = NULL, method="nlminb", arnoldWied=FALSE, 
 control = list(), pars, scaleU=FALSE, verbose=NULL, legacy=FALSE,
 se.lambda=TRUE, returnHcov=FALSE, pWOrder=250, tol.Hcov=1.0e-10)
\method{summary}{gmsar}(object, correlation = FALSE, Hausman=FALSE, ...)
GMargminImage(obj, lambdaseq, s2seq)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a symbolic description of the model to be fit. The details 
of model specification are given for \code{lm()}}
  \item{data}{an optional data frame containing the variables in the model. 
By default the variables are taken from the environment which the function 
is called.}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{na.action}{a function (default \code{na.fail}), can also be \code{na.omit} or \code{na.exclude} with consequences for residuals and fitted values - in these cases the weights list will be subsetted to remove NAs in the data. It may be necessary to set zero.policy to TRUE because this subsetting may create no-neighbour observations. Note that only weights lists created without using the glist argument to \code{nb2listw} may be subsetted.}
  \item{zero.policy}{default NULL, use global option value; if TRUE assign zero to the lagged value of zones without 
neighbours, if FALSE (default) assign NA - causing \code{GMerrorsar()} to terminate with an error}
  \item{method}{default \code{"nlminb"}, or optionally a method passed to \code{optim} to use an alternative optimizer}
  \item{arnoldWied}{default FALSE}
  \item{control}{A list of control parameters. See details in \code{\link[stats]{optim}} or \code{\link[stats]{nlminb}}.}
  \item{pars}{starting values for \eqn{\lambda}{lambda} and \eqn{\sigma^2}{sigma squared} for GMM optimisation, if missing (default), approximated from initial OLS model as the autocorrelation coefficient corrected for weights style and model sigma squared}
  \item{scaleU}{Default FALSE: scale the OLS residuals before computing the moment matrices; only used if the \code{pars} argument is missing}
  \item{verbose}{default NULL, use global option value; if TRUE, reports function values during optimization.}
  \item{legacy}{default FALSE - compute using the unfiltered values of the response and right hand side variables; if TRUE - compute the fitted value and residuals from the spatially filtered model using the spatial error parameter}
  \item{se.lambda}{default TRUE, use the analytical method described in \url{http://econweb.umd.edu/~prucha/STATPROG/OLS/desols.pdf}}
  \item{returnHcov}{default FALSE, return the Vo matrix for a spatial Hausman test}
  \item{tol.Hcov}{the tolerance for computing the Vo matrix (default=1.0e-10)}
  \item{pWOrder}{default 250, if returnHcov=TRUE, pass this order to \code{powerWeights} as the power series maximum limit}
  \item{object, obj}{\code{gmsar} object from \code{GMerrorsar}}
  \item{correlation}{logical; (default=FALSE), TRUE not available}
  \item{Hausman}{if TRUE, the results of the Hausman test for error models are reported}
  \item{\dots}{\code{summary} arguments passed through}
  \item{lambdaseq}{if given, an increasing sequence of lambda values for gridding}
  \item{s2seq}{if given, an increasing sequence of sigma squared values for gridding}
}
\details{
  When the control list is set with care, the function will converge to values close to the ML estimator without requiring computation of the Jacobian, the most resource-intensive part of ML estimation. 

Note that the fitted() function for the output object assumes that the response 
variable may be reconstructed as the sum of the trend, the signal, and the
noise (residuals). Since the values of the response variable are known,
their spatial lags are used to calculate signal components (Cressie 1993, p. 564). This differs from other software, including GeoDa, which does not use knowledge of the response 
variable in making predictions for the fitting data.

The \code{GMargminImage} may be used to visualize the shape of the surface of the argmin function used to find lambda.
}
\value{
  A list object of class \code{gmsar}
  \item{type}{"ERROR"}
  \item{lambda}{simultaneous autoregressive error coefficient}
  \item{coefficients}{GMM coefficient estimates}
  \item{rest.se}{GMM coefficient standard errors}
  \item{s2}{GMM residual variance}
  \item{SSE}{sum of squared GMM errors}
  \item{parameters}{number of parameters estimated}
  \item{lm.model}{the \code{lm} object returned when estimating for \eqn{\lambda=0}{lambda=0}}
  \item{call}{the call used to create this object}
  \item{residuals}{GMM residuals}
  \item{lm.target}{the \code{lm} object returned for the GMM fit}
  \item{fitted.values}{Difference between residuals and response variable}
  \item{formula}{model formula}
  \item{aliased}{if not NULL, details of aliased variables}
  \item{zero.policy}{zero.policy for this model}
  \item{vv}{list of internal bigG and litg components for testing optimisation surface}
  \item{optres}{object returned by optimizer}
  \item{pars}{start parameter values for optimisation}
  \item{Hcov}{Spatial DGP covariance matrix for Hausman test if available}
  \item{legacy}{input choice of unfiltered or filtered values}
  \item{lambda.se}{value computed if input argument TRUE}
  \item{arnoldWied}{were Arnold-Wied moments used}
  \item{GMs2}{GM argmin sigma squared}
  \item{scaleU}{input choice of scaled OLS residuals}
  \item{vcov}{variance-covariance matrix of regression coefficients}
  \item{na.action}{(possibly) named vector of excluded or omitted observations if non-default na.action argument used}
}
\references{Kelejian, H. H., and Prucha, I. R., 1999. A Generalized Moments Estimator for the Autoregressive Parameter in a Spatial Model. International Economic Review, 40, pp. 509--533; Cressie, N. A. C. 1993 \emph{Statistics for spatial data}, Wiley, New York.

Roger Bivand, Gianfranco Piras (2015). Comparing Implementations of Estimation Methods for Spatial Econometrics. \emph{Journal of Statistical Software}, 63(18), 1-36. \url{https://www.jstatsoft.org/v63/i18/}.
}
\author{Luc Anselin and Roger Bivand}

\seealso{\code{\link[stats]{optim}}, \code{\link[stats]{nlminb}}, \code{\link{errorsarlm}}}
\examples{
data(oldcol)
COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="W"), method="eigen")
summary(COL.errW.eig, Hausman=TRUE)
COL.errW.GM <- GMerrorsar(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="W"), returnHcov=TRUE)
summary(COL.errW.GM, Hausman=TRUE)
aa <- GMargminImage(COL.errW.GM)
levs <- quantile(aa$z, seq(0, 1, 1/12))
image(aa, breaks=levs, xlab="lambda", ylab="s2")
points(COL.errW.GM$lambda, COL.errW.GM$s2, pch=3, lwd=2)
contour(aa, levels=signif(levs, 4), add=TRUE)
COL.errW.GM1 <- GMerrorsar(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="W"))
summary(COL.errW.GM1)
if (require(foreign, quietly=TRUE)) {
example(NY_data, package="spData")
esar1f <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, family="SAR", method="eigen")
summary(esar1f)
esar1gm <- GMerrorsar(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
 data=nydata, listw=listw_NY)
summary(esar1gm)
esar1gm1 <- GMerrorsar(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
 data=nydata, listw=listw_NY, method="Nelder-Mead")
summary(esar1gm1)
}
}
\keyword{spatial}

