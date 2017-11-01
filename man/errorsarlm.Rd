% Copyright 2002-16 by Roger S. Bivand
\name{errorsarlm}
\alias{errorsarlm}
\alias{lmSLX}
\alias{create_WX}
%\alias{sar.error.f}
%\alias{sar.error.f.s}
\title{Spatial simultaneous autoregressive error model estimation}
\description{
  Maximum likelihood estimation of spatial simultaneous autoregressive
error models of the form:

\deqn{y = X \beta + u, u = \lambda W u + \varepsilon}{y = X beta + u, u = lambda W u + e}

where \eqn{\lambda}{lambda} is found by \code{optimize()} first, and \eqn{\beta}{beta} and other parameters by generalized least squares subsequently. With one of the sparse matrix methods, larger numbers of observations can be handled, but the \code{interval=} argument may need be set when the weights are not row-standardised. When \code{etype} is \dQuote{emixed}, a so-called spatial Durbin error model is fitted, while \code{lmSLX} fits an \code{lm} model augmented with the spatially lagged RHS variables, including the lagged intercept when the spatial weights are not row-standardised. \code{create_WX} creates spatially lagged RHS variables, and is exposed for use in model fitting functions.
}
\usage{
errorsarlm(formula, data=list(), listw, na.action, weights=NULL,
 etype="error", method="eigen", quiet=NULL, zero.policy=NULL,
 interval = NULL, tol.solve=1.0e-10, trs=NULL, control=list())
lmSLX(formula, data = list(), listw, na.action, weights=NULL, zero.policy=NULL)
create_WX(x, listw, zero.policy=NULL, prefix="")
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{formula}{a symbolic description of the model to be fit. The details 
of model specification are given for \code{lm()}}
  \item{data}{an optional data frame containing the variables in the model. 
By default the variables are taken from the environment which the function 
is called.}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{na.action}{a function (default \code{options("na.action")}), can also be \code{na.omit} or \code{na.exclude} with consequences for residuals and fitted values - in these cases the weights list will be subsetted to remove NAs in the data. It may be necessary to set zero.policy to TRUE because this subsetting may create no-neighbour observations. Note that only weights lists created without using the glist argument to \code{nb2listw} may be subsetted.}
  \item{weights}{an optional vector of weights to be used in the fitting process. Non-NULL weights can be used to indicate that different observations have different variances (with the values in weights being inversely proportional to the variances); or equivalently, when the elements of weights are positive integers w_i, that each response y_i is the mean of w_i unit-weight observations (including the case that there are w_i observations equal to y_i and the data have been summarized) - \code{\link{lm}}}
  \item{etype}{default "error", may be set to "emixed" to include the spatially lagged independent variables added to X; when "emixed", the lagged intercept is dropped for spatial weights style "W", that is row-standardised weights, but otherwise included}
  \item{method}{"eigen" (default) - the Jacobian is computed as the product 
of (1 - rho*eigenvalue) using \code{eigenw}, and "spam" or "Matrix_J" for strictly symmetric weights lists of styles "B" and "C", or made symmetric by similarity (Ord, 1975, Appendix C) if possible for styles "W" and "S", using code from the spam package or Matrix package to calculate the determinant; \dQuote{Matrix} and \dQuote{spam_update} provide updating Cholesky decomposition methods; "LU" provides an alternative sparse matrix decomposition approach. In addition, there are "Chebyshev" and Monte Carlo "MC" approximate log-determinant methods; the Smirnov/Anselin (2009) trace approximation is available as "moments". Three methods: "SE_classic", "SE_whichMin", and "SE_interp" are provided experimentally, the first to attempt to emulate the behaviour of Spatial Econometrics toolbox ML fitting functions. All use grids of log determinant values, and the latter two attempt to ameliorate some features of "SE_classic".
}
  \item{quiet}{default NULL, use !verbose global option value; if FALSE, reports function values during optimization.}
  \item{zero.policy}{default NULL, use global option value; if TRUE assign zero to the lagged value of zones without 
neighbours, if FALSE assign NA - causing \code{errorsarlm()} to terminate with an error}
  \item{interval}{default is NULL, search interval for autoregressive parameter}
  \item{tol.solve}{the tolerance for detecting linear dependencies in the columns of matrices to be inverted - passed to \code{solve()} (default=1.0e-10). This may be used if necessary to extract coefficient standard errors (for instance lowering to 1e-12), but errors in \code{solve()} may constitute indications of poorly scaled variables: if the variables have scales differing much from the autoregressive coefficient, the values in this matrix may be very different in scale, and inverting such a matrix is analytically possible by definition, but numerically unstable; rescaling the RHS variables alleviates this better than setting tol.solve to a very small value}
  \item{trs}{default NULL, if given, a vector of powered spatial weights matrix traces output by \code{trW}; when given, insert the asymptotic analytical values into the numerical Hessian instead of the approximated values; may be used to get around some problems raised when the numerical Hessian is poorly conditioned, generating NaNs in subsequent operations. When using the numerical Hessian to get the standard error of lambda, it is very strongly advised that trs be given, as the parts of fdHess corresponding to the regression coefficients are badly approximated, affecting the standard error of lambda; the coefficient correlation matrix is unusable}
 \item{control}{list of extra control arguments - see section below}
 \item{x}{model matrix to be lagged}
 \item{prefix}{default empty string, may be \dQuote{lag} in some cases}
}
\details{
The asymptotic standard error of \eqn{\lambda}{lambda} is only computed when
method=eigen, because the full matrix operations involved would be costly
for large n typically associated with the choice of method="spam" or
"Matrix".  The same applies to the coefficient covariance matrix. Taken
as the asymptotic matrix from the literature, it is typically badly
scaled, being block-diagonal, and with the elements involving \eqn{\lambda}{lambda}
being very small, while other parts of the matrix can be very large
(often many orders of magnitude in difference). It often happens that
the \code{tol.solve} argument needs to be set to a smaller value than
the default, or the RHS variables can be centred or reduced in range.

Note that the fitted() function for the output object assumes that the response 
variable may be reconstructed as the sum of the trend, the signal, and the
noise (residuals). Since the values of the response variable are known,
their spatial lags are used to calculate signal components (Cressie 1993, p. 564). This differs from other software, including GeoDa, which does not use knowledge of the response 
variable in making predictions for the fitting data.
}

\section{Control arguments}{
\describe{
  \item{tol.opt:}{the desired accuracy of the optimization - passed to \code{optimize()} (default=square root of double precision machine tolerance, a larger root may be used needed, see help(boston) for an example)}
  \item{returnHcov:}{default TRUE, return the Vo matrix for a spatial Hausman test}
  \item{pWOrder:}{default 250, if returnHcov=TRUE and the method is not \dQuote{eigen}, pass this order to \code{powerWeights} as the power series maximum limit}
  \item{fdHess:}{default NULL, then set to (method != "eigen") internally; use \code{fdHess} to compute an approximate Hessian using finite differences when using sparse matrix methods; used to make a coefficient covariance matrix when the number of observations is large; may be turned off to save resources if need be}
  \item{optimHess:}{default FALSE, use \code{fdHess} from \pkg{nlme}, if TRUE, use \code{optim} to calculate Hessian at optimum}
  \item{optimHessMethod:}{default \dQuote{optimHess}, may be \dQuote{nlm} or one of the \code{optim} methods}
  \item{LAPACK:}{default FALSE; logical value passed to \code{qr} in the SSE log likelihood function}
  \item{compiled_sse:}{default FALSE; logical value used in the log likelihood function to choose compiled code for computing SSE}
  \item{Imult:}{default 2; used for preparing the Cholesky decompositions for updating in the Jacobian function}
  \item{super:}{if NULL (default), set to FALSE to use a simplicial decomposition for the sparse Cholesky decomposition and method \dQuote{Matrix_J}, set to  \code{as.logical(NA)} for method \dQuote{Matrix}, if TRUE, use a supernodal decomposition}
  \item{cheb_q:}{default 5; highest power of the approximating polynomial for the Chebyshev approximation}
  \item{MC_p:}{default 16; number of random variates}
  \item{MC_m:}{default 30; number of products of random variates matrix and spatial weights matrix}
  \item{spamPivot:}{default \dQuote{MMD}, alternative \dQuote{RCM}}
  \item{in_coef}{default 0.1, coefficient value for initial Cholesky decomposition in \dQuote{spam_update}}
  \item{type}{default \dQuote{MC}, used with method \dQuote{moments}; alternatives \dQuote{mult} and \dQuote{moments}, for use if \code{trs} is missing, \code{\link{trW}}}
  \item{correct}{default TRUE, used with method \dQuote{moments} to compute the Smirnov/Anselin correction term}
  \item{trunc}{default TRUE, used with method \dQuote{moments} to truncate the Smirnov/Anselin correction term}
  \item{SE_method}{default \dQuote{LU}, may be \dQuote{MC}}
  \item{nrho}{default 200, as in SE toolbox; the size of the first stage lndet grid; it may be reduced to for example 40}
  \item{interpn}{default 2000, as in SE toolbox; the size of the second stage lndet grid}
  \item{small_asy}{default TRUE; if the method is not \dQuote{eigen}, use asymmetric covariances rather than numerical Hessian ones if n <= small}
  \item{small}{default 1500; threshold number of observations for asymmetric covariances when the method is not \dQuote{eigen}}
  \item{SElndet}{default NULL, may be used to pass a pre-computed SE toolbox style matrix of coefficients and their lndet values to the "SE_classic" and "SE_whichMin" methods}
  \item{LU_order}{default FALSE; used in \dQuote{LU_prepermutate}, note warnings given for \code{lu} method}
  \item{pre_eig}{default NULL; may be used to pass a pre-computed vector of eigenvalues}
}}

\value{
  A list object of class \code{sarlm}
  \item{type}{"error"}
  \item{lambda}{simultaneous autoregressive error coefficient}
  \item{coefficients}{GLS coefficient estimates}
  \item{rest.se}{GLS coefficient standard errors (are equal to asymptotic
standard errors)}
  \item{LL}{log likelihood value at computed optimum}
  \item{s2}{GLS residual variance}
  \item{SSE}{sum of squared GLS errors}
  \item{parameters}{number of parameters estimated}
  \item{logLik_lm.model}{Log likelihood of the linear model for \eqn{\lambda=0}{lambda=0}}
  \item{AIC_lm.model}{AIC of the linear model for \eqn{\lambda=0}{lambda=0}}
%  \item{lm.model}{the \code{lm} object returned when estimating for \eqn{\lambda=0}{lambda=0}}
  \item{coef_lm.model}{coefficients of the linear model for \eqn{\lambda=0}{lambda=0}}
  \item{tarX}{model matrix of the GLS model}
  \item{tary}{response of the GLS model}
  \item{y}{response of the linear model for \eqn{\lambda=0}{lambda=0}}
  \item{X}{model matrix of the linear model for \eqn{\lambda=0}{lambda=0}}
  \item{method}{the method used to calculate the Jacobian}
  \item{call}{the call used to create this object}
  \item{residuals}{GLS residuals}
%  \item{lm.target}{the \code{lm} object returned for the GLS fit}
  \item{opt}{object returned from numerical optimisation}
  \item{fitted.values}{Difference between residuals and response variable}
  \item{ase}{TRUE if method=eigen}
%  \item{formula}{model formula}
  \item{se.fit}{Not used yet}
  \item{lambda.se}{if ase=TRUE, the asymptotic standard error of \eqn{\lambda}{lambda}}
  \item{LMtest}{NULL for this model}
  \item{aliased}{if not NULL, details of aliased variables}
  \item{LLNullLlm}{Log-likelihood of the null linear model}
  \item{Hcov}{Spatial DGP covariance matrix for Hausman test if available}
  \item{interval}{line search interval}
  \item{fdHess}{finite difference Hessian}
  \item{optimHess}{\code{optim} or \code{fdHess} used}
  \item{insert}{logical; is TRUE, asymptotic values inserted in fdHess where feasible}
  \item{timings}{processing timings}
  \item{f_calls}{number of calls to the log likelihood function during optimization}
  \item{hf_calls}{number of calls to the log likelihood function during numerical Hessian computation}
  \item{intern_classic}{a data frame of detval matrix row choices used by the SE toolbox classic method}
  \item{zero.policy}{zero.policy for this model}
  \item{na.action}{(possibly) named vector of excluded or omitted observations if non-default na.action argument used}
  \item{weights}{weights used in model fitting}
  \item{emixedImps}{for \dQuote{emixed} models, a  list of three impact matrixes (impacts and standard errors) for direct, indirect and total impacts; total impacts calculated using gmodels::estimable}

The internal sar.error.* functions return the value of the log likelihood function at \eqn{\lambda}{lambda}.

The \code{lmSLX} function returns an \dQuote{lm} object with a \dQuote{mixedImps} list of three impact matrixes (impacts and standard errors) for direct, indirect and total impacts; total impacts calculated using gmodels::estimable.
}
\references{Cliff, A. D., Ord, J. K. 1981 \emph{Spatial processes}, Pion;
Ord, J. K. 1975 Estimation methods for models of spatial interaction,
\emph{Journal of the American Statistical Association}, 70, 120-126;
Anselin, L. 1988 \emph{Spatial econometrics: methods and models.}
(Dordrecht: Kluwer); Anselin, L. 1995 SpaceStat, a software program for
the analysis of spatial data, version 1.80. Regional Research Institute,
West Virginia University, Morgantown, WV;
Anselin L, Bera AK (1998) Spatial dependence in linear regression models
with an introduction to spatial econometrics. In: Ullah A, Giles DEA
(eds) Handbook of applied economic statistics. Marcel Dekker, New York,
pp. 237-289; Cressie, N. A. C. 1993 \emph{Statistics for spatial data}, Wiley, New York; LeSage J and RK Pace (2009) Introduction to Spatial Econometrics. CRC Press, Boca Raton.

Roger Bivand, Gianfranco Piras (2015). Comparing Implementations of Estimation Methods for Spatial Econometrics. \emph{Journal of Statistical Software}, 63(18), 1-36. \url{https://www.jstatsoft.org/v63/i18/}.

Bivand, R. S., Hauke, J., and Kossowski, T. (2013). Computing the Jacobian in Gaussian spatial autoregressive models: An illustrated comparison of available methods. \emph{Geographical Analysis}, 45(2), 150-179.
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{lm}}, \code{\link{lagsarlm}}, \code{\link{similar.listw}}, \code{\link{summary.sarlm}}, \code{\link{predict.sarlm}},
\code{\link{residuals.sarlm}}, \code{\link{do_ldet}}, \code{\link[gmodels]{estimable}}
}

\examples{
data(oldcol)
lw <- nb2listw(COL.nb, style="W")
COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="eigen", quiet=FALSE)
summary(COL.errW.eig, correlation=TRUE)
ev <- eigenw(similar.listw(lw))
COL.errW.eig_ev <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="eigen", control=list(pre_eig=ev))
all.equal(coefficients(COL.errW.eig), coefficients(COL.errW.eig_ev))
COL.errB.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="B"), method="eigen", quiet=FALSE)
summary(COL.errB.eig, correlation=TRUE)
W <- as(nb2listw(COL.nb), "CsparseMatrix")
trMatc <- trW(W, type="mult")
COL.errW.M <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="Matrix", quiet=FALSE, trs=trMatc)
summary(COL.errW.M, correlation=TRUE)
COL.SDEM.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="eigen", etype="emixed")
summary(COL.SDEM.eig, correlation=TRUE)
summary(impacts(COL.SDEM.eig))
summary(impacts(COL.SDEM.eig), adjust_k=TRUE)
COL.SLX <- lmSLX(CRIME ~ INC + HOVAL, data=COL.OLD, listw=lw)
summary(COL.SLX)
summary(impacts(COL.SLX))
COL.SLX <- lmSLX(CRIME ~ INC + HOVAL + I(HOVAL^2), data=COL.OLD, listw=lw)
summary(COL.SLX)
COL.SLX <- lmSLX(CRIME ~ INC, data=COL.OLD, listw=lw)
\donttest{
crds <- cbind(COL.OLD$X, COL.OLD$Y)
mdist <- sqrt(sum(diff(apply(crds, 2, range))^2))
dnb <- dnearneigh(crds, 0, mdist)
dists <- nbdists(dnb, crds)
f <- function(x, form, data, dnb, dists, verbose) {
  glst <- lapply(dists, function(d) 1/(d^x))
  lw <- nb2listw(dnb, glist=glst, style="B")
  res <- logLik(lmSLX(form=form, data=data, listw=lw))
  if (verbose) cat("power:", x, "logLik:", res, "\n")
  res
}
opt <- optimize(f, interval=c(0.1, 4), form=CRIME ~ INC + HOVAL,
 data=COL.OLD, dnb=dnb, dists=dists, verbose=TRUE, maximum=TRUE)
glst <- lapply(dists, function(d) 1/(d^opt$maximum))
lw <- nb2listw(dnb, glist=glst, style="B")
SLX <- lmSLX(CRIME ~ INC + HOVAL, data=COL.OLD, listw=lw)
summary(SLX)
summary(impacts(SLX))
}
NA.COL.OLD <- COL.OLD
NA.COL.OLD$CRIME[20:25] <- NA
COL.err.NA <- errorsarlm(CRIME ~ INC + HOVAL, data=NA.COL.OLD,
 nb2listw(COL.nb), na.action=na.exclude)
COL.err.NA$na.action
COL.err.NA
resid(COL.err.NA)
\donttest{
lw <- nb2listw(COL.nb, style="W")
system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="eigen"))
ocoef <- coefficients(COL.errW.eig)
system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="eigen", control=list(LAPACK=FALSE)))
all.equal(ocoef, coefficients(COL.errW.eig))
system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="eigen", control=list(compiled_sse=TRUE)))
all.equal(ocoef, coefficients(COL.errW.eig))
system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="Matrix_J", control=list(super=TRUE)))
all.equal(ocoef, coefficients(COL.errW.eig))
system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="Matrix_J", control=list(super=FALSE)))
all.equal(ocoef, coefficients(COL.errW.eig))
system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="Matrix_J", control=list(super=as.logical(NA))))
all.equal(ocoef, coefficients(COL.errW.eig))
system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="Matrix", control=list(super=TRUE)))
all.equal(ocoef, coefficients(COL.errW.eig))
system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="Matrix", control=list(super=FALSE)))
all.equal(ocoef, coefficients(COL.errW.eig))
system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="Matrix", control=list(super=as.logical(NA))))
all.equal(ocoef, coefficients(COL.errW.eig))
system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="spam", control=list(spamPivot="MMD")))
all.equal(ocoef, coefficients(COL.errW.eig))
system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="spam", control=list(spamPivot="RCM")))
all.equal(ocoef, coefficients(COL.errW.eig))
system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="spam_update", control=list(spamPivot="MMD")))
all.equal(ocoef, coefficients(COL.errW.eig))
system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="spam_update", control=list(spamPivot="RCM")))
all.equal(ocoef, coefficients(COL.errW.eig))
}
}
\keyword{spatial}
