\name{spdep-deprecated}
\alias{spdep-deprecated}
\alias{lextrB}
\alias{lextrW}
\alias{lextrS}
\alias{l_max}
\alias{griffith_sone}
\alias{subgraph_eigenw}
\alias{mom_calc}
\alias{mom_calc_int2}
\alias{stsls}
\alias{print.stsls}
\alias{print.summary.stsls}
\alias{summary.stsls}
\alias{residuals.stsls}
\alias{coef.stsls}
\alias{deviance.stsls}
\alias{impacts.stsls}
\alias{GMerrorsar}
\alias{residuals.gmsar}
\alias{deviance.gmsar}
\alias{coef.gmsar}
\alias{fitted.gmsar}
\alias{print.gmsar}
\alias{summary.gmsar}
\alias{print.summary.gmsar}
\alias{Hausman.test.gmsar}
\alias{GMargminImage}
\alias{gstsls}
\alias{impacts.lagmess}
\alias{lagmess}
\alias{print.lagmess}
\alias{print.summary.lagmess}
\alias{summary.lagmess}
\alias{residuals.lagmess}
\alias{deviance.lagmess}
\alias{coef.lagmess}
\alias{fitted.lagmess}
\alias{logLik.lagmess}
\alias{ME}
\alias{print.ME_res}
\alias{fitted.ME_res}
\alias{SpatialFiltering}
\alias{print.SFResult}
\alias{fitted.SFResult}
\alias{LR.sarlm}
\alias{LR1.sarlm}
\alias{Wald1.sarlm}
\alias{Hausman.test}
\alias{Hausman.test.sarlm}
\alias{logLik.sarlm}
\alias{as.spam.listw}
\alias{as_dgRMatrix_listw}
\alias{as_dsTMatrix_listw}
\alias{as_dsCMatrix_I}
\alias{as_dsCMatrix_IrW}
\alias{Jacobian_W}
\alias{coerce,listw,CsparseMatrix-method}
\alias{coerce,listw,RsparseMatrix-method}
\alias{coerce,listw,symmetricMatrix-method}
\alias{powerWeights}
\alias{impacts.lagmess}
\alias{plot.lagImpact}
\alias{print.lagImpact}
\alias{summary.lagImpact}
\alias{print.summary.lagImpact}
\alias{HPDinterval.lagImpact}
\alias{can.be.simmed}
\alias{eigenw}
\alias{similar.listw}
\alias{do_ldet}
\alias{jacobianSetup}
\alias{eigen_setup}
\alias{eigen_pre_setup}
\alias{mcdet_setup}
\alias{cheb_setup}
\alias{spam_setup}
\alias{spam_update_setup}
\alias{Matrix_setup}
\alias{Matrix_J_setup}
\alias{LU_setup}
\alias{LU_prepermutate_setup}
\alias{moments_setup}
\alias{SE_classic_setup}
\alias{SE_whichMin_setup}
\alias{SE_interp_setup}
\alias{MCMCsamp}
\alias{MCMCsamp.spautolm}
\alias{MCMCsamp.sarlm}
\alias{spautolm}
\alias{residuals.spautolm}
\alias{deviance.spautolm}
\alias{coef.spautolm}
\alias{fitted.spautolm}
\alias{print.spautolm}
\alias{summary.spautolm}
\alias{LR1.spautolm}
\alias{logLik.spautolm}
\alias{print.summary.spautolm}
\alias{spBreg_sac}
\alias{impacts.MCMC_sar_g}
\alias{impacts.MCMC_sem_g}
\alias{impacts.MCMC_sac_g}
\alias{spBreg_err}
\alias{spBreg_lag}
\alias{predict.SLX}
\alias{lmSLX}
\alias{impacts.SLX}
\alias{create_WX}
\alias{anova.sarlm}
\alias{bptest.sarlm}
\alias{errorsarlm}
\alias{impacts}
\alias{impacts.sarlm}
\alias{intImpacts}
\alias{lagsarlm}
\alias{predict.sarlm}
\alias{print.sarlm.pred}
\alias{as.data.frame.sarlm.pred}
\alias{residuals.sarlm}
\alias{deviance.sarlm}
\alias{coef.sarlm}
\alias{vcov.sarlm}
\alias{fitted.sarlm}
\alias{sacsarlm}
\alias{summary.sarlm}
\alias{print.sarlm}
\alias{print.summary.sarlm}


\title{Deprecated Functions in Package \pkg{spdep}}

\description{
  These functions are provided for compatibility with older versions of
  \pkg{spdep} only, and may be defunct as soon as the next release. The
  functions have been moved to the \pkg{spatialreg} package.
}

\usage{
lextrB(lw, zero.policy = TRUE, control = list())
lextrW(lw, zero.policy=TRUE, control=list())
lextrS(lw, zero.policy=TRUE, control=list())
l_max(lw, zero.policy=TRUE, control=list())
griffith_sone(P, Q, type="rook")
subgraph_eigenw(nb, glist=NULL, style="W", zero.policy=NULL, quiet=NULL)
mom_calc(lw, m)
mom_calc_int2(is, m, nb, weights, Card)
stsls(formula, data = list(), listw, zero.policy = NULL,
 na.action = na.fail, robust = FALSE, HC=NULL, legacy=FALSE, W2X = TRUE)
\method{impacts}{stsls}(obj, \dots, tr, R = NULL, listw = NULL, evalues=NULL,
 tol = 1e-06, empirical = FALSE, Q=NULL)
GMerrorsar(formula, data = list(), listw, na.action = na.fail,
 zero.policy = NULL, method="nlminb", arnoldWied=FALSE, 
 control = list(), pars=NULL, scaleU=FALSE, verbose=NULL, legacy=FALSE,
 se.lambda=TRUE, returnHcov=FALSE, pWOrder=250, tol.Hcov=1.0e-10)
\method{summary}{gmsar}(object, correlation = FALSE, Hausman=FALSE, ...)
GMargminImage(obj, lambdaseq, s2seq)
gstsls(formula, data = list(), listw, listw2 = NULL, na.action = na.fail, 
    zero.policy = NULL, pars=NULL, scaleU=FALSE, control = list(), 
    verbose=NULL, method="nlminb", robust=FALSE, legacy=FALSE, W2X=TRUE) 
\method{impacts}{gmsar}(obj, \dots, n = NULL, tr = NULL, R = NULL,
 listw = NULL, evalues=NULL, tol = 1e-06, empirical = FALSE, Q=NULL)
\method{Hausman.test}{gmsar}(object, ..., tol=NULL)
lagmess(formula, data = list(), listw, zero.policy = NULL, na.action = na.fail,
 q = 10, start = -2.5, control=list(), method="BFGS", verbose=NULL,
 use_expm=FALSE)
ME(formula, data=list(), family = gaussian, weights, offset,
 na.action=na.fail, listw=NULL, alpha=0.05, nsim=99, verbose=NULL,
 stdev=FALSE, zero.policy = NULL)
SpatialFiltering(formula, lagformula=NULL, data=list(), na.action=na.fail,
 nb=NULL, glist = NULL, style = "C", zero.policy = NULL, tol = 0.1,
 zerovalue = 1e-04, ExactEV = FALSE, symmetric = TRUE, alpha=NULL,
 alternative="two.sided", verbose=NULL)
LR.sarlm(x, y)
\method{logLik}{sarlm}(object, ...)
LR1.sarlm(object)
Wald1.sarlm(object)
\method{Hausman.test}{sarlm}(object, ..., tol=NULL)
as.spam.listw(listw)
as_dgRMatrix_listw(listw)
as_dsTMatrix_listw(listw)
as_dsCMatrix_I(n)
as_dsCMatrix_IrW(W, rho)
Jacobian_W(W, rho)
powerWeights(W, rho, order=250, X, tol=.Machine$double.eps^(3/5))
\method{plot}{lagImpact}(x, ..., choice="direct", trace=FALSE, density=TRUE)
\method{print}{lagImpact}(x, ..., reportQ=NULL)
\method{summary}{lagImpact}(object, ..., zstats=FALSE, short=FALSE, reportQ=NULL)
\method{HPDinterval}{lagImpact}(obj, prob = 0.95, ..., choice="direct")
intImpacts(rho, beta, P, n, mu, Sigma, irho, drop2beta, bnames, interval,
 type, tr, R, listw, evalues, tol, empirical, Q, icept, iicept, p, mess=FALSE,
 samples=NULL, zero_fill = NULL, dvars = NULL)
can.be.simmed(listw)
eigenw(listw, quiet=NULL)
similar.listw(listw)
do_ldet(coef, env, which=1)
jacobianSetup(method, env, con, pre_eig=NULL, trs=NULL, interval=NULL, which=1)
cheb_setup(env, q=5, which=1)
mcdet_setup(env, p=16, m=30, which=1)
eigen_setup(env, which=1)
eigen_pre_setup(env, pre_eig, which=1)
spam_setup(env, pivot="MMD", which=1)
spam_update_setup(env, in_coef=0.1, pivot="MMD", which=1)
Matrix_setup(env, Imult, super=as.logical(NA), which=1)
Matrix_J_setup(env, super=FALSE, which=1)
LU_setup(env, which=1)
LU_prepermutate_setup(env, coef=0.1, order=FALSE, which=1)
moments_setup(env, trs=NULL, m, p, type="MC", correct=TRUE, trunc=TRUE, eq7=TRUE, which=1)
SE_classic_setup(env, SE_method="LU", p=16, m=30, nrho=200, interpn=2000,
 interval=c(-1,0.999), SElndet=NULL, which=1)
SE_whichMin_setup(env, SE_method="LU", p=16, m=30, nrho=200, interpn=2000,
 interval=c(-1,0.999), SElndet=NULL, which=1)
SE_interp_setup(env, SE_method="LU", p=16, m=30, nrho=200,
 interval=c(-1,0.999), which=1)
MCMCsamp(object, mcmc = 1L, verbose = NULL, ...)
\method{MCMCsamp}{spautolm}(object, mcmc = 1L, verbose = NULL, ...,
 burnin = 0L, scale=1, listw, control = list())
\method{MCMCsamp}{sarlm}(object, mcmc = 1L, verbose = NULL, ...,
    burnin=0L, scale=1, listw, listw2=NULL, control=list())
spautolm(formula, data = list(), listw, weights,
 na.action, family = "SAR", method="eigen", verbose = NULL, trs=NULL,
 interval=NULL, zero.policy = NULL, tol.solve=.Machine$double.eps,
 llprof=NULL, control=list())
\method{summary}{spautolm}(object, correlation = FALSE, adj.se=FALSE,
 Nagelkerke=FALSE, ...)
spBreg_sac(formula, data = list(), listw, listw2=NULL, na.action, 
    Durbin, type, zero.policy=NULL, control=list())
\method{impacts}{MCMC_sar_g}(obj, ..., tr=NULL, listw=NULL, evalues=NULL, Q=NULL)
\method{impacts}{MCMC_sem_g}(obj, ..., tr=NULL, listw=NULL, evalues=NULL, Q=NULL)
\method{impacts}{MCMC_sac_g}(obj, ..., tr=NULL, listw=NULL, evalues=NULL, Q=NULL)
spBreg_err(formula, data = list(), listw, na.action, Durbin, etype,
    zero.policy=NULL, control=list())
spBreg_lag(formula, data = list(), listw, na.action, Durbin, type,
    zero.policy=NULL, control=list())
\method{predict}{SLX}(object, newdata, listw, zero.policy=NULL, ...)
lmSLX(formula, data = list(), listw, na.action, weights=NULL,
 Durbin=TRUE, zero.policy=NULL)
\method{impacts}{SLX}(obj, ...)
create_WX(x, listw, zero.policy=NULL, prefix="")
\method{anova}{sarlm}(object, ...)
bptest.sarlm(object, varformula=NULL, studentize = TRUE, data=list())
errorsarlm(formula, data=list(), listw, na.action, weights=NULL,
 Durbin, etype, method="eigen", quiet=NULL, zero.policy=NULL,
 interval = NULL, tol.solve=1.0e-10, trs=NULL, control=list())
\method{impacts}{sarlm}(obj, \dots, tr, R = NULL, listw = NULL, evalues=NULL,
 useHESS = NULL, tol = 1e-06, empirical = FALSE, Q=NULL)
lagsarlm(formula, data = list(), listw, 
	na.action, Durbin, type, method="eigen", quiet=NULL, 
	zero.policy=NULL, interval=NULL, tol.solve=1.0e-10, trs=NULL, 
	control=list())
\method{predict}{sarlm}(object, newdata = NULL, listw = NULL, pred.type = "TS", all.data = FALSE,
 zero.policy = NULL, legacy = TRUE, legacy.mixed = FALSE, power = NULL, order = 250,
 tol = .Machine$double.eps^(3/5), spChk = NULL, ...)
\method{print}{sarlm.pred}(x, ...)
\method{as.data.frame}{sarlm.pred}(x, ...)
\method{residuals}{sarlm}(object, ...)
\method{deviance}{sarlm}(object, ...)
\method{coef}{sarlm}(object, ...)
\method{vcov}{sarlm}(object, ...)
\method{fitted}{sarlm}(object, ...)
sacsarlm(formula, data = list(), listw, listw2 = NULL, na.action, Durbin, type,
 method = "eigen", quiet = NULL, zero.policy = NULL, tol.solve = 1e-10,
 llprof=NULL, interval1=NULL, interval2=NULL, trs1=NULL, trs2=NULL,
 control = list())
\method{summary}{sarlm}(object, correlation = FALSE, Nagelkerke = FALSE, Hausman=FALSE, adj.se=FALSE, ...)
\method{print}{sarlm}(x, ...)
\method{print}{summary.sarlm}(x, digits = max(5, .Options$digits - 3),
	signif.stars = FALSE, ...)

}

\arguments{
  \item{lw}{a binary symmetric \code{listw} object from, for example, \code{nb2listw} with style \dQuote{B} for \code{lextrB}, style \dQuote{W} for \code{lextrW} and style \dQuote{S} for \code{lextrS}; for \code{l_max}, the object may be asymmetric and does not have to be binary}
  \item{zero.policy}{default NULL, use global option value; if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA}
  \item{control}{a list of control arguments}
  \item{quiet}{default NULL, use global !verbose option value; set to FALSE for short summary}
  \item{P}{number of columns in the grid (number of units in a horizontal axis direction)}
  \item{Q}{number of rows in the grid (number of units in a vertical axis direction.)}
  \item{type}{\dQuote{rook} or \dQuote{queen}}
  \item{nb}{an object of class \code{nb}}
  \item{glist}{list of general weights corresponding to neighbours}
  \item{style}{\code{style} can take values \dQuote{W}, \dQuote{B}, \dQuote{C}, \dQuote{U}, \dQuote{minmax} and \dQuote{S}}
  \item{m}{The number of powers; must be an even number for \sQuote{type}=\dQuote{moments} (default changed from 100 to 30 (2010-11-17))}
  \item{is}{(used internally only in \code{mom_calc_int2} for \sQuote{type}=\dQuote{moments} on a cluster)}
%  \item{nb}{(used internally only in \code{mom_calc_int2} for \sQuote{type}=\dQuote{moments} on a cluster)}
  \item{weights}{(used internally only in \code{mom_calc_int2} for \sQuote{type}=\dQuote{moments} on a cluster)}
  \item{Card}{(used internally only in \code{mom_calc_int2} for \sQuote{type}=\dQuote{moments} on a cluster)}
  \item{formula}{a symbolic description of the model to be fit. The details 
of model specification are given for \code{lm()}}
  \item{data}{an optional data frame containing the variables in the model. 
By default the variables are taken from the environment which the function 
is called.}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{na.action}{a function (default \code{na.fail}), can also be \code{na.omit} or \code{na.exclude} with consequences for residuals and fitted values - in these cases the weights list will be subsetted to remove NAs in the data. It may be necessary to set zero.policy to TRUE because this subsetting may create no-neighbour observations. Note that only weights lists created without using the glist argument to \code{nb2listw} may be subsetted.}
  \item{robust}{default FALSE, if TRUE, apply a heteroskedasticity correction to the coefficients covariances}
  \item{HC}{default NULL, if \code{robust} is TRUE, assigned \dQuote{HC0}, may take values \dQuote{HC0} or \dQuote{HC1} for White estimates or MacKinnon-White estimates respectively}
 \item{legacy}{the argument chooses between two implementations of the robustness correction: default FALSE - use the estimate of Omega only in the White consistent estimator of the variance-covariance matrix, if TRUE, use the original implementation which runs a GLS using the estimate of Omega, and yields different coefficient estimates as well - see example below}
 \item{W2X}{default TRUE, if FALSE only WX are used as instruments in the spatial two stage least squares; until release 0.4-60, only WX were used - see example below }
  \item{obj}{A spatial regression object created by \code{lagsarlm}, \code{lagmess} or by \code{lmSLX}; in \code{HPDinterval.lagImpact}, a lagImpact object}
  \item{...}{Arguments passed through to methods in the \pkg{coda} package}
  \item{tr}{A vector of traces of powers of the spatial weights matrix created using \code{trW}, for approximate impact measures; if not given, \code{listw} must be given for exact measures (for small to moderate spatial weights matrices); the traces must be for the same spatial weights as were used in fitting the spatial regression, and must be row-standardised}
  \item{evalues}{vector of eigenvalues of spatial weights matrix for impacts calculations}
  \item{n}{defaults to \code{length(obj$residuals)}; in the method for \code{gmsar} objects it may be used in panel settings to compute the impacts for cross-sectional weights only, suggested by Angela Parenti}
  \item{R}{If given, simulations are used to compute distributions for the impact measures, returned as \code{mcmc} objects; the objects are used for convenience but are not output by an MCMC process}
%  \item{useHESS}{Use the Hessian approximation (if available) even if the asymptotic coefficient covariance matrix is available; used for comparing methods}
  \item{tol}{Argument passed to \code{mvrnorm}: tolerance (relative to largest variance) for numerical lack of positive-definiteness in the coefficient covariance matrix}
  \item{empirical}{Argument passed to \code{mvrnorm} (default FALSE): if true, the coefficients and their covariance matrix specify the empirical not population mean and covariance matrix}
  \item{listw2}{a \code{listw} object created for example by \code{nb2listw}, if not given, set to the same spatial weights as the listw argument}
  \item{pars}{starting values for \eqn{\lambda}{lambda} and \eqn{\sigma^2}{sigma squared} for GMM optimisation, 
  if missing (default), approximated from initial 2sls model as the autocorrelation coefficient corrected for weights style 
  and model sigma squared}
  \item{scaleU}{Default FALSE: scale the OLS residuals before computing the moment matrices; only used if the \code{pars} argument is missing}
  \item{method}{default \code{"nlminb"}, or optionally a method passed to \code{optim} to use an alternative optimizer}
  \item{arnoldWied}{default FALSE}
  \item{returnHcov}{default FALSE, return the Vo matrix for a spatial Hausman test}
  \item{tol.Hcov}{the tolerance for computing the Vo matrix (default=1.0e-10)}
  \item{pWOrder}{default 250, if returnHcov=TRUE, pass this order to \code{powerWeights} as the power series maximum limit}
  \item{lambdaseq}{if given, an increasing sequence of lambda values for gridding}
  \item{s2seq}{if given, an increasing sequence of sigma squared values for gridding}
  \item{object}{\code{gmsar} object from \code{GMerrorsar}}
  \item{correlation}{logical; (default=FALSE), TRUE not available}
  \item{Hausman}{if TRUE, the results of the Hausman test for error models are reported}
  \item{se.lambda}{default TRUE, use the analytical method described in \url{http://econweb.umd.edu/~prucha/STATPROG/OLS/desols.pdf}}
  \item{verbose}{default NULL, use global option value; if TRUE, reports function values during optimization.}
  \item{q}{default 10; number of powers of the spatial weights to use}
  \item{start}{starting value for numerical optimization, should be a small negative number}
  \item{use_expm}{default FALSE; if TRUE use \code{expm::expAtv} instead of a truncated power series of W}
  \item{family}{a description of the error distribution and link function to
          be used in the model}
%  \item{weights}{an optional vector of weights to be used in the fitting process}
  \item{offset}{this can be used to specify an a priori known component to be included in the linear predictor during fitting}
  \item{alpha}{used as a stopping rule to choose all eigenvectors up to and including the one with a p-value exceeding alpha}
  \item{nsim}{number of permutations for permutation bootstrap for finding p-values}
  \item{stdev}{if TRUE, p-value calculated from bootstrap permutation standard deviate using \code{pnorm} with alternative="greater", if FALSE the Hope-type p-value}
  \item{lagformula}{An extra one-sided formula to be used when a spatial lag representation is desired; the intercept is excluded within the function if present because it is part of the formula argument, but excluding it explicitly in the lagformula argument in the presence of factors generates a collinear model matrix}
  \item{zerovalue}{eigenvectors with eigenvalues of an absolute value smaller than zerovalue will be excluded in eigenvector search}
  \item{ExactEV}{Set ExactEV=TRUE to use exact expectations and variances rather than the expectation and variance of Moran's I from the previous iteration, default FALSE}
  \item{symmetric}{Should the spatial weights matrix be forced to symmetry, default TRUE}
  \item{alternative}{a character string specifying the alternative hypothesis, must be one of greater, less or two.sided (default).}
  \item{x}{a \code{logLik} object or an object for which a \code{logLik()} function exists}
  \item{y}{a \code{logLik} object or an object for which a \code{logLik()} function exists}
%  \item{object}{a \code{sarlm} object from \code{lagsarlm()} or \code{errorsarlm()}}
%  \item{\dots}{further arguments passed to or from other methods}
%  \item{tol}{\code{tol} argument passed to \code{solve}, default NULL}
  \item{W}{a \code{dsTMatrix} object created using \code{as_dsTMatrix_listw} from a symmetric \code{listw} object}
  \item{rho}{spatial regression coefficient}
  \item{order}{Power series maximum limit}
  \item{X}{A numerical matrix}
%  \item{tol}{Tolerance for convergence of power series}
  \item{choice}{One of three impacts: direct, indirect, or total}
  \item{trace}{Argument passed to \code{plot.mcmc}: plot trace plots}
  \item{density}{Argument passed to \code{plot.mcmc}: plot density plots}
  \item{prob}{Argument passed to \code{HPDinterval.mcmc}: a numeric scalar in the interval (0,1) giving the target probability content of the intervals}
%  \item{adjust_k}{default TRUE if SDEM else FALSE, adjust internal OLS SDEM standard errors by dividing by n rather than (n-k) (default changed and bug fixed after 0.7-8; standard errors now ML in SDEM summary and impacts summary and identical - for SLX use FALSE)}
  \item{
beta, 
mu, Sigma, irho, drop2beta, bnames,
    interval, 
 icept, iicept, p, mess, samples, zero_fill, dvars}{internal arguments shared inside impacts methods}
  \item{reportQ}{default NULL; if TRUE and \code{Q} given as an argument to \code{impacts}, report impact components}
%  \item{x, object}{lagImpact objects created by \code{impacts} methods}
  \item{zstats}{default FALSE, if TRUE, also return z-values and p-values for the impacts based on the simulations}
  \item{short}{default FALSE, if TRUE passed to the print summary method to omit printing of the mcmc summaries}
  \item{coef}{spatial coefficient value}
  \item{env}{environment containing pre-computed objects, fixed after assignment in setup functions}
  \item{which}{default 1; if 2, use second listw object}
%  \item{method}{string value, used by \code{jacobianSetup} to choose method}
  \item{con}{control list passed from model fitting function and parsed in \code{jacobianSetup} to set environment variables for method-specific setup}
  \item{pre_eig}{pre-computed eigenvalues of length n}
%  \item{q}{Chebyshev approximation order; default in calling spdep functions is 5, here it cannot be missing and does not have a default}
%  \item{p}{Monte Carlo approximation number of random normal variables; default calling spdep functions is 16, here it cannot be missing and does not have a default}
%  \item{m}{Monte Carlo approximation number of series terms; default in calling spdep functions is 30, here it cannot be missing and does not have a default; \code{m} serves the same purpose in the moments method}
  \item{pivot}{default \dQuote{MMD}, may also be \dQuote{RCM} for Cholesky decompisition using spam}
  \item{in_coef}{fill-in initiation coefficient value, default 0.1}
  \item{Imult}{see \code{\link[Matrix]{Cholesky}}; numeric scalar which defaults to zero. The matrix that is decomposed is A+m*I where m is the value of Imult and I is the identity matrix of order ncol(A). Default in calling spdep functions is 2, here it cannot be missing and does not have a default, but is rescaled for binary weights matrices in proportion to the maximim row sum in those calling functions}
  \item{super}{see \code{\link[Matrix]{Cholesky}}; logical scalar indicating is a supernodal decomposition should be created.  The alternative is a simplicial decomposition. Default in calling spdep functions is FALSE for \dQuote{Matrix_J} and \code{as.logical(NA)} for \dQuote{Matrix}.  Setting it to NA leaves the choice to a CHOLMOD-internal heuristic}
%  \item{order}{default FALSE; used in LU_prepermutate, note warnings given for \code{lu} method}
  \item{trs, trs1, trs2}{A numeric vector of \code{m} traces, as from \code{trW}}
%  \item{type}{moments trace type, see \code{\link{trW}}}
  \item{correct}{default TRUE: use Smirnov correction term, see \code{\link{trW}}}
  \item{trunc}{default TRUE: truncate Smirnov correction term, see \code{\link{trW}}}
  \item{eq7}{default TRUE}{use equation 7 in Smirnov and Anselin (2009), if FALSE no unit root correction}
  \item{SE_method}{default \dQuote{LU}, alternatively \dQuote{MC}; underlying lndet method to use for generating SE toolbox emulation grid}
  \item{nrho}{default 200, number of lndet values in first stage SE toolbox emulation grid}
  \item{interval1, interval2}{default c(-1,0.999) if interval argument NULL, bounds for SE toolbox emulation grid}
  \item{interpn}{default 2000, number of lndet values to interpolate in second stage SE toolbox emulation grid}
  \item{SElndet}{default NULL, used to pass a pre-computed two-column matrix of coefficient values and corresponding interpolated lndet values}
  \item{mcmc}{The number of MCMC iterations after burnin}
%  \item{verbose}{default NULL, use global option value; if TRUE, reports progress}
%  \item{\dots}{Arguments passed through}
  \item{burnin}{The number of burn-in iterations for the sampler}
  \item{scale}{a positive scale parameter}
%  \item{listw, listw2}{\code{listw} objects created for example by \code{nb2listw}; should be the same object(s) used for fitting the model}
%  \item{control}{list of extra control arguments - see \code{\link{spautolm}}}
%  \item{formula}{a symbolic description of the model to be fit. The details 
%of model specification are given for \code{lm()}}
%  \item{data}{an optional data frame containing the variables in the model. 
%By default the variables are taken from the environment which the function 
%is called.}
%  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
%  \item{weights}{an optional vector of weights to be used in the fitting process}
%  \item{na.action}{a function (default \code{options("na.action")}), can also be \code{na.omit} or \code{na.exclude} with consequences for residuals and fitted values - in these cases the weights list will be subsetted to remove NAs in the data. Note that only weights lists created without using the glist argument to \code{nb2listw} may be subsetted.}
%  \item{family}{character string: either \code{"SAR"} or \code{"CAR"} for simultaneous or conditional autoregressions; \code{"SMA"} for spatial moving average added thanks to Jielai Ma - \code{"SMA"} is only implemented for method=\code{"eigen"} because it necessarily involves dense matrices}
%  \item{method}{character string: default \code{"eigen"} for use of dense matrices, \code{"Matrix_J"} for sparse matrices (restricted to spatial weights symmetric or similar to symmetric) using methods in the Matrix package; \dQuote{Matrix} provides updating Cholesky decomposition methods. Values of method may also include "LU", which provides an alternative sparse matrix decomposition approach, and the "Chebyshev" and Monte Carlo "MC" approximate log-determinant methods.}
%  \item{verbose}{default NULL, use global option value; if TRUE, reports function values during optimization.}
%  \item{trs}{default NULL, if given, a vector of powered spatial weights matrix traces output by \code{trW}; when given, used in some Jacobian methods}
%  \item{interval}{search interval for autoregressive parameter when not using method="eigen"; default is c(-1,0.999), \code{optimize} will reset NA/NaN to a bound and gives a warning when the interval is poorly set; method="Matrix" will attempt to search for an appropriate interval, if find\_interval=TRUE (fails on some platforms)}
%  \item{zero.policy}{default NULL, use global option value; Include list of no-neighbour observations in output if TRUE --- otherwise zero.policy is handled within the listw argument}
  \item{tol.solve}{the tolerance for detecting linear dependencies in the columns of matrices to be inverted - passed to \code{solve()} (default=double precision machine tolerance). Errors in \code{solve()} may constitute indications of poorly scaled variables: if the variables have scales differing much from the autoregressive coefficient, the values in this matrix may be very different in scale, and inverting such a matrix is analytically possible by definition, but numerically unstable; rescaling the RHS variables alleviates this better than setting tol.solve to a very small value}
  \item{llprof}{default NULL, can either be an integer, to divide the feasible range into llprof points, or a sequence of spatial coefficient values, at which to evaluate the likelihood function}
% \item{control}{list of extra control arguments - see section below}
%  \item{object}{\code{spautolm} object from \code{spautolm}}
%  \item{correlation}{logical; if 'TRUE', the correlation matrix of the estimated parameters is returned and printed (default=FALSE)}
  \item{adj.se}{if TRUE, adjust the coefficient standard errors for the number of fitted coefficients}
  \item{Nagelkerke}{if TRUE, the Nagelkerke pseudo R-squared is reported}
  \item{Durbin}{default FALSE (spatial lag model); if TRUE, full spatial Durbin model; if a formula object, the subset of explanatory variables to lag}
  \item{etype}{(use the \sQuote{Durbin=} argument - retained for backwards compatibility only) default "error", may be set to "emixed" to include the spatially lagged independent variables added to X; when "emixed", the lagged intercept is dropped for spatial weights style "W", that is row-standardised weights, but otherwise included}
  \item{newdata}{data frame in which to predict --- if NULL, predictions are
for the data on which the model was fitted. Should have row names corresponding to region.id. If row names are exactly the same than the ones used for training, it uses in-sample predictors for forecast.  See \sQuote{Details}}
 \item{prefix}{default empty string, may be \dQuote{lag} in some cases}
  \item{varformula}{a formula describing only the potential explanatory variables for the variance (no dependent variable needed). By default the same explanatory variables are taken as in the main regression model}
  \item{studentize}{logical. If set to \code{TRUE} Koenker's studentized
   version of the test statistic will be used.}
  \item{useHESS, pred.type, all.data, legacy.mixed, power, spChk, digits, signif.stars}{other arguments in deprecated functions}
}

\details{
 Model-fitting functions and functions supporting model fitting are being moved to the \pkg{spatialreg} package.
}


\seealso{
  \code{\link[base]{Deprecated}}
}

\examples{
\dontrun{
data(boston, package="spData")
ab.listb <- nb2listw(boston.soi, style="B")
er <- range(eigenw(ab.listb))
er
res_1 <- lextrB(ab.listb)
c(res_1)
run <- FALSE
if (require("RSpectra", quietly=TRUE)) run <- TRUE
if (run) {
B <- as(ab.listb, "CsparseMatrix")
eigs(B, k=1, which="SR")$values
}
if (run) {
eigs(B, k=1, which="LR")$values
}
k5 <- knn2nb(knearneigh(boston.utm, k=5))
c(l_max(nb2listw(k5, style="B")))
max(Re(eigenw(nb2listw(k5, style="B"))))
c(l_max(nb2listw(k5, style="C")))
max(Re(eigenw(nb2listw(k5, style="C"))))
ab.listw <- nb2listw(boston.soi, style="W")
er <- range(eigenw(similar.listw(ab.listw)))
er
res_1 <- lextrW(ab.listw)
c(res_1)
if (run) {
B <- as(similar.listw(ab.listw), "CsparseMatrix")
eigs(B, k=1, which="SR")$values
}
if (run) {
eigs(B, k=1, which="LR")$values
}
ab.listw <- nb2listw(boston.soi, style="S")
er <- range(eigenw(similar.listw(ab.listw)))
er
res_1 <- lextrS(ab.listw)
c(res_1)
if (run) {
B <- as(similar.listw(ab.listw), "CsparseMatrix")
eigs(B, k=1, which="SR")$values
}
if (run) {
eigs(B, k=1, which="LR")$values
}
rg <- cell2nb(ncol=7, nrow=7, type="rook")
rg_eig <- eigenw(nb2listw(rg, style="B"))
rg_GS <- griffith_sone(P=7, Q=7, type="rook")
all.equal(rg_eig, rg_GS)
# subgraphs
data(oldcol)
crds <- cbind(COL.OLD$X, COL.OLD$Y)
k3 <- knn2nb(knearneigh(crds, k=3))
nc <- n.comp.nb(k3)
nc$nc
table(nc$comp.id)
k3eig <- eigenw(nb2listw(k3, style="W"))
k3eigSG <- subgraph_eigenw(k3, style="W")
all.equal(sort(k3eig), k3eigSG)

data(oldcol)
COL.lag.eig <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, nb2listw(COL.nb))
summary(COL.lag.eig, correlation=TRUE)
COL.lag.stsls <- stsls(CRIME ~ INC + HOVAL, data=COL.OLD, nb2listw(COL.nb))
summary(COL.lag.stsls, correlation=TRUE)
COL.lag.stslsW <- stsls(CRIME ~ INC + HOVAL, data=COL.OLD, nb2listw(COL.nb), W2X=FALSE)
summary(COL.lag.stslsW, correlation=TRUE)
COL.lag.stslsR <- stsls(CRIME ~ INC + HOVAL, data=COL.OLD, nb2listw(COL.nb),
robust=TRUE, W2X=FALSE)
summary(COL.lag.stslsR, correlation=TRUE)
COL.lag.stslsRl <- stsls(CRIME ~ INC + HOVAL, data=COL.OLD, nb2listw(COL.nb),
robust=TRUE, legacy=TRUE, W2X=FALSE)
summary(COL.lag.stslsRl, correlation=TRUE)
data(boston, package="spData")
gp2a <- stsls(log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + I(RM^2) +
  AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT),
 data=boston.c, nb2listw(boston.soi))
summary(gp2a)
columbus <- st_read(system.file("shapes/columbus.shp", package="spData")[1], quiet=TRUE)
col.gal.nb <- read.gal(system.file("weights/columbus.gal", package="spData")[1])
listw <- nb2listw(col.gal.nb)
ev <- eigenw(listw)
W <- as(listw, "CsparseMatrix")
trMatc <- trW(W, type="mult")
lobj1 <- stsls(CRIME ~ INC + HOVAL, columbus, listw)
loobj1 <- impacts(lobj1, R=200, tr=trMatc)
summary(loobj1, zstats=TRUE, short=TRUE)
loobj2 <- impacts(lobj1, R=200, evalues=ev)
summary(loobj2, zstats=TRUE, short=TRUE)
library(coda)
HPDinterval(loobj1)
lobj1r <- stsls(CRIME ~ INC + HOVAL, columbus, listw, robust=TRUE)
loobj1r <- impacts(lobj1r, tr=trMatc, R=200)
summary(loobj1r, zstats=TRUE, short=TRUE)
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
nydata <- st_read(system.file("shapes/NY8_bna_utm18.gpkg", package="spData")[1], quiet=TRUE)
suppressMessages(nyadjmat <- as.matrix(foreign::read.dbf(system.file(
 "misc/nyadjwts.dbf", package="spData")[1])[-1]))
suppressMessages(ID <- as.character(names(foreign::read.dbf(system.file(
 "misc/nyadjwts.dbf", package="spData")[1]))[-1]))
identical(substring(ID, 2, 10), substring(as.character(nydata$AREAKEY), 2, 10))
nyadjlw <- mat2listw(nyadjmat, as.character(nydata$AREAKEY))
listw_NY <- nb2listw(nyadjlw$neighbours, style="B")
esar1f <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, family="SAR", method="eigen")
summary(esar1f)
esar1gm <- GMerrorsar(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
 data=nydata, listw=listw_NY)
summary(esar1gm)
esar1gm1 <- GMerrorsar(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
 data=nydata, listw=listw_NY, method="Nelder-Mead")
summary(esar1gm1)
data(baltimore, package="spData")
baltimore$AGE <- ifelse(baltimore$AGE < 1, 1, baltimore$AGE)
lw <- nb2listw(knn2nb(knearneigh(cbind(baltimore$X, baltimore$Y), k=7)))
obj1 <- lm(log(PRICE) ~ PATIO + log(AGE) + log(SQFT),
 data=baltimore)
lm.morantest(obj1, lw)
lm.LMtests(obj1, lw, test="all")
system.time(obj2 <- lagmess(log(PRICE) ~ PATIO + log(AGE) + log(SQFT), data=baltimore, listw=lw))
summary(obj2)
system.time(obj2a <- lagmess(log(PRICE) ~ PATIO + log(AGE) + log(SQFT), data=baltimore, listw=lw,
 use_expm=TRUE))
summary(obj2a)
obj3 <- lagsarlm(log(PRICE) ~ PATIO + log(AGE) + log(SQFT), data=baltimore, listw=lw)
summary(obj3)
data(boston, package="spData")
lw <- nb2listw(boston.soi)
gp2 <- lagsarlm(log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + I(RM^2)
 +  AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT),
 data=boston.c, lw, method="Matrix")
summary(gp2)
gp2a <- lagmess(CMEDV ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + I(RM^2)
 +  AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT),
 data=boston.c, lw)
summary(gp2a)
columbus <- st_read(system.file("shapes/columbus.shp", package="spData")[1], quiet=TRUE)
col.gal.nb <- read.gal(system.file("weights/columbus.gal", package="spData")[1])
lmbase <- lm(CRIME ~ INC + HOVAL, data=columbus)
lagcol <- SpatialFiltering(CRIME ~ 1, ~ INC + HOVAL, data=columbus,
 nb=col.gal.nb, style="W", alpha=0.1, verbose=TRUE)
lagcol
lmlag <- lm(CRIME ~ INC + HOVAL + fitted(lagcol), data=columbus)
anova(lmbase, lmlag)
set.seed(123)
lagcol1 <- ME(CRIME ~ INC + HOVAL, data=columbus, family="gaussian",
 listw=nb2listw(col.gal.nb), alpha=0.1, verbose=TRUE)
lagcol1
lmlag1 <- lm(CRIME ~ INC + HOVAL + fitted(lagcol1), data=columbus)
anova(lmbase, lmlag1)
set.seed(123)
lagcol2 <- ME(CRIME ~ INC + HOVAL, data=columbus, family="gaussian",
 listw=nb2listw(col.gal.nb), alpha=0.1, stdev=TRUE, verbose=TRUE)
lagcol2
lmlag2 <- lm(CRIME ~ INC + HOVAL + fitted(lagcol2), data=columbus)
anova(lmbase, lmlag2)
NA.columbus <- columbus
NA.columbus$CRIME[20:25] <- NA
COL.ME.NA <- ME(CRIME ~ INC + HOVAL, data=NA.columbus, family="gaussian",
 listw=nb2listw(col.gal.nb), alpha=0.1, stdev=TRUE, verbose=TRUE,
 na.action=na.exclude)
COL.ME.NA$na.action
summary(lm(CRIME ~ INC + HOVAL + fitted(COL.ME.NA), data=NA.columbus,
 na.action=na.exclude))
#nc.sids <- st_read(system.file("shapes/sids.shp", package="spData")[1], quiet=TRUE)
#rn <- as.character(nc.sids$FIPS)
#ncCC89_nb <- read.gal(system.file("weights/ncCC89.gal", package="spData")[1],
# region.id=rn)
#ncCR85_nb <- read.gal(system.file("weights/ncCR85.gal", package="spData")[1],
# region.id=rn)
#glmbase <- glm(SID74 ~ 1, data=nc.sids, offset=log(BIR74),
# family="poisson")
#set.seed(123)
#MEpois1 <- ME(SID74 ~ 1, data=nc.sids, offset=log(BIR74),
# family="poisson", listw=nb2listw(ncCR85_nb, style="B"), alpha=0.2, verbose=TRUE)
#MEpois1
#glmME <- glm(SID74 ~ 1 + fitted(MEpois1), data=nc.sids, offset=log(BIR74),
# family="poisson")
#anova(glmME, test="Chisq")
#anova(glmbase, glmME, test="Chisq")
data(hopkins, package="spData")
hopkins_part <- hopkins[21:36,36:21]
hopkins_part[which(hopkins_part > 0, arr.ind=TRUE)] <- 1
hopkins.rook.nb <- cell2nb(16, 16, type="rook")
glmbase <- glm(c(hopkins_part) ~ 1, family="binomial")
set.seed(123)
MEbinom1 <- ME(c(hopkins_part) ~ 1, family="binomial",
 listw=nb2listw(hopkins.rook.nb, style="B"), alpha=0.2, verbose=TRUE)
glmME <- glm(c(hopkins_part) ~ 1 + fitted(MEbinom1), family="binomial")
anova(glmME, test="Chisq")
anova(glmbase, glmME, test="Chisq")
columbus <- st_read(system.file("shapes/columbus.shp", package="spData")[1], quiet=TRUE)
col.gal.nb <- read.gal(system.file("weights/columbus.gal", package="spData")[1])
lmbase <- lm(CRIME ~ INC + HOVAL, data=columbus)
sarcol <- SpatialFiltering(CRIME ~ INC + HOVAL, data=columbus,
 nb=col.gal.nb, style="W", ExactEV=TRUE)
sarcol
lmsar <- lm(CRIME ~ INC + HOVAL + fitted(sarcol), data=columbus)
lmsar
anova(lmbase, lmsar)
lm.morantest(lmsar, nb2listw(col.gal.nb))
lagcol <- SpatialFiltering(CRIME ~ 1, ~ INC + HOVAL - 1, data=columbus,
 nb=col.gal.nb, style="W")
lagcol
lmlag <- lm(CRIME ~ INC + HOVAL + fitted(lagcol), data=columbus)
lmlag
anova(lmbase, lmlag)
lm.morantest(lmlag, nb2listw(col.gal.nb))
NA.columbus <- columbus
NA.columbus$CRIME[20:25] <- NA
COL.SF.NA <- SpatialFiltering(CRIME ~ INC + HOVAL, data=NA.columbus,
 nb=col.gal.nb, style="W", na.action=na.exclude)
COL.SF.NA$na.action
summary(lm(CRIME ~ INC + HOVAL + fitted(COL.SF.NA), data=NA.columbus,
 na.action=na.exclude))
columbus <- st_read(system.file("shapes/columbus.shp", package="spData")[1], quiet=TRUE)
col.gal.nb <- read.gal(system.file("weights/columbus.gal", package="spData")[1])
mixed <- lagsarlm(CRIME ~ HOVAL + INC, data=columbus, nb2listw(col.gal.nb),
  type="mixed")
error <- errorsarlm(CRIME ~ HOVAL + INC, data=columbus, nb2listw(col.gal.nb))
LR.sarlm(mixed, error)
Hausman.test(error)
columbus <- st_read(system.file("shapes/columbus.shp", package="spData")[1], quiet=TRUE)
col.gal.nb <- read.gal(system.file("weights/columbus.gal", package="spData")[1])
col.listw <- nb2listw(col.gal.nb)
if (require("spam", quietly=TRUE)) {
col.sp <- as.spam.listw(col.listw)
str(col.sp)
}
suppressMessages(nyadjmat <- as.matrix(foreign::read.dbf(system.file(
 "misc/nyadjwts.dbf", package="spData")[1])[-1]))
nyadjlw <- mat2listw(nyadjmat)
listw_NY <- nb2listw(nyadjlw$neighbours, style="B")
W_C <- as(listw_NY, "CsparseMatrix")
W_R <- as(listw_NY, "RsparseMatrix")
W_S <- as(listw_NY, "symmetricMatrix")
n <- nrow(W_S)
I <- Diagonal(n)
rho <- 0.1
c(determinant(I - rho * W_S, logarithm=TRUE)$modulus)
sum(log(1 - rho * eigenw(listw_NY)))
nW <- - W_S
nChol <- Cholesky(nW, Imult=8)
n * log(rho) + (2 * c(determinant(update(nChol, nW, 1/rho))$modulus))
nb7rt <- cell2nb(7, 7, torus=TRUE)
x <- matrix(sample(rnorm(500*length(nb7rt))), nrow=length(nb7rt))
W <- as(nb2listw(nb7rt), "CsparseMatrix")
system.time(ee <- powerWeights(W, rho=0.9, X=x))
#all.equal(e, as(ee, "matrix"), check.attributes=FALSE)
system.time(e <- invIrM(nb7rt, rho=0.98, method="solve", feasible=NULL) \%*\% x)
system.time(ee <- powerWeights(W, rho=0.98, X=x))
str(attr(ee, "internal"))
all.equal(e, as(ee, "matrix"), check.attributes=FALSE)
system.time(ee <- powerWeights(W, rho=0.98, order=1000, X=x))
all.equal(e, as(ee, "matrix"), check.attributes=FALSE)
nb60rt <- cell2nb(60, 60, torus=TRUE)
W <- as(nb2listw(nb60rt), "CsparseMatrix")
set.seed(1)
x <- matrix(rnorm(dim(W)[1]), ncol=1)
system.time(ee <- powerWeights(W, rho=0.3, X=x))
str(as(ee, "matrix"))
obj <- errorsarlm(as(ee, "matrix")[,1] ~ 1, listw=nb2listw(nb60rt), method="Matrix")
coefficients(obj)
data(oldcol)
W.eig <- eigenw(nb2listw(COL.nb, style="W"))
1/range(W.eig)
S.eig <- eigenw(nb2listw(COL.nb, style="S"))
1/range(S.eig)
B.eig <- eigenw(nb2listw(COL.nb, style="B"))
1/range(B.eig)
# cases for intrinsically asymmetric weights
crds <- cbind(COL.OLD$X, COL.OLD$Y)
k3 <- knn2nb(knearneigh(crds, k=3))
is.symmetric.nb(k3)
k3eig <- eigenw(nb2listw(k3, style="W"))
is.complex(k3eig)
rho <- 0.5
Jc <- sum(log(1 - rho * k3eig))
# complex eigenvalue Jacobian
Jc
W <- as(nb2listw(k3, style="W"), "CsparseMatrix")
I <- diag(length(k3))
Jl <- sum(log(abs(diag(slot(lu(I - rho * W), "U")))))
# LU Jacobian equals complex eigenvalue Jacobian
Jl
all.equal(Re(Jc), Jl)
# wrong value if only real part used
Jr <- sum(log(1 - rho * Re(k3eig)))
Jr
all.equal(Jr, Jl)
# construction of Jacobian from complex conjugate pairs (Jan Hauke)
Rev <- Re(k3eig)[which(Im(k3eig) == 0)]
# real eigenvalues
Cev <- k3eig[which(Im(k3eig) != 0)]
pCev <- Cev[Im(Cev) > 0]
# separate complex conjugate pairs
RpCev <- Re(pCev)
IpCev <- Im(pCev)
# reassemble Jacobian
Jc1 <- sum(log(1 - rho*Rev)) + sum(log((1 - rho * RpCev)^2 + (rho^2)*(IpCev^2)))
all.equal(Re(Jc), Jc1)
# impact of omitted complex part term in real part only Jacobian
Jc2 <- sum(log(1 - rho*Rev)) + sum(log((1 - rho * RpCev)^2))
all.equal(Jr, Jc2)
# trace of asymmetric (WW) and crossprod of complex eigenvalues for APLE
sum(diag(W \%*\% W))
crossprod(k3eig)
# analytical regular grid eigenvalues
run <- FALSE
if (require("RSpectra", quietly=TRUE)) run <- TRUE
if (run) {
rg <- cell2nb(ncol=7, nrow=7, type="rook")
B <- as(nb2listw(rg, style="B"), "CsparseMatrix")
res1 <- eigs(B, k=1, which="LR")$values
resn <- eigs(B, k=1, which="SR")$values
print(Re(c(resn, res1)))
}
if (run) {
rg_eig <- eigenw(nb2listw(rg, style="B"))
print(all.equal(range(Re(rg_eig)), c(resn, res1))) 
}
if (run) {
lw <- nb2listw(rg, style="W")
rg_eig <- eigenw(similar.listw(lw))
print(range(Re(rg_eig)))
}
if (run) {
W  <- as(lw, "CsparseMatrix")
print(Re(c(eigs(W, k=1, which="SR")$values, eigs(W, k=1, which="LR")$values)))
}
data(oldcol)
COL.W <- nb2listw(COL.nb, style="W")
COL.S <- nb2listw(COL.nb, style="S")
sum(log(1 - 0.5 * eigenw(COL.W)))
sum(log(1 - 0.5 * eigenw(similar.listw(COL.W))))
W_J <- as(as_dsTMatrix_listw(similar.listw(COL.W)), "CsparseMatrix")
I <- as_dsCMatrix_I(dim(W_J)[1])
c(determinant(I - 0.5 * W_J, logarithm=TRUE)$modulus)
sum(log(1 - 0.5 * eigenw(COL.S)))
sum(log(1 - 0.5 * eigenw(similar.listw(COL.S))))
W_J <- as(as_dsTMatrix_listw(similar.listw(COL.S)), "CsparseMatrix")
c(determinant(I - 0.5 * W_J, logarithm=TRUE)$modulus)
data(boston, package="spData")
lw <- nb2listw(boston.soi)
can.sim <- spdep:::can.be.simmed(lw)
env <- new.env(parent=globalenv())
assign("listw", lw, envir=env)
assign("can.sim", can.sim, envir=env)
assign("similar", FALSE, envir=env)
assign("verbose", FALSE, envir=env)
assign("family", "SAR", envir=env)
eigen_setup(env)
get("similar", envir=env)
do_ldet(0.5, env)
rm(env)
env <- new.env(parent=globalenv())
assign("listw", lw, envir=env)
assign("can.sim", can.sim, envir=env)
assign("similar", FALSE, envir=env)
assign("verbose", FALSE, envir=env)
assign("family", "SAR", envir=env)
assign("n", length(boston.soi), envir=env)
eigen_pre_setup(env, pre_eig=eigenw(similar.listw(lw)))
do_ldet(0.5, env)
rm(env)
env <- new.env(parent=globalenv())
assign("listw", lw, envir=env)
assign("can.sim", can.sim, envir=env)
assign("similar", FALSE, envir=env)
assign("family", "SAR", envir=env)
assign("n", length(boston.soi), envir=env)
Matrix_setup(env, Imult=2, super=FALSE)
get("similar", envir=env)
do_ldet(0.5, env)
rm(env)
env <- new.env(parent=globalenv())
assign("listw", lw, envir=env)
assign("n", length(boston.soi), envir=env)
assign("can.sim", can.sim, envir=env)
assign("similar", FALSE, envir=env)
assign("family", "SAR", envir=env)
spam_setup(env)
get("similar", envir=env)
do_ldet(0.5, env)
rm(env)
env <- new.env(parent=globalenv())
assign("listw", lw, envir=env)
assign("n", length(boston.soi), envir=env)
assign("similar", FALSE, envir=env)
assign("family", "SAR", envir=env)
LU_setup(env)
get("similar", envir=env)
do_ldet(0.5, env)
rm(env)
env <- new.env(parent=globalenv())
assign("listw", lw, envir=env)
assign("n", length(boston.soi), envir=env)
assign("similar", FALSE, envir=env)
assign("family", "SAR", envir=env)
LU_prepermutate_setup(env)
get("similar", envir=env)
do_ldet(0.5, env)
rm(env)
env <- new.env(parent=globalenv())
assign("listw", lw, envir=env)
assign("similar", FALSE, envir=env)
assign("family", "SAR", envir=env)
cheb_setup(env, q=5)
get("similar", envir=env)
do_ldet(0.5, env)
rm(env)
env <- new.env(parent=globalenv())
assign("listw", lw, envir=env)
assign("n", length(boston.soi), envir=env)
assign("similar", FALSE, envir=env)
assign("family", "SAR", envir=env)
set.seed(12345)
mcdet_setup(env, p=16, m=30)
get("similar", envir=env)
do_ldet(0.5, env)
rm(env)
nydata <- st_read(system.file("shapes/NY8_bna_utm18.gpkg", package="spData")[1], quiet=TRUE)
suppressMessages(nyadjmat <- as.matrix(foreign::read.dbf(system.file(
 "misc/nyadjwts.dbf", package="spData")[1])[-1]))
suppressMessages(ID <- as.character(names(foreign::read.dbf(system.file(
 "misc/nyadjwts.dbf", package="spData")[1]))[-1]))
identical(substring(ID, 2, 10), substring(as.character(nydata$AREAKEY), 2, 10))
nyadjlw <- mat2listw(nyadjmat)
listw_NY <- nb2listw(nyadjlw$neighbours, style="B")
esar1f <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, family="SAR", method="eigen")
summary(esar1f)
res <- MCMCsamp(esar1f, mcmc=5000, burnin=500, listw=listw_NY)
summary(res)
#esar1fw <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
# listw=listw_NY, weights=POP8, family="SAR", method="eigen")
#summary(esar1fw)
#res <- MCMCsamp(esar1fw, mcmc=5000, burnin=500, listw=listw_NY)
#summary(res)
ecar1f <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, family="CAR", method="eigen")
summary(ecar1f)
res <- MCMCsamp(ecar1f, mcmc=5000, burnin=500, listw=listw_NY)
summary(res)
#esar1fw <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
# listw=listw_NY, weights=POP8, family="SAR", method="eigen")
#summary(esar1fw)
#res <- MCMCsamp(esar1fw, mcmc=5000, burnin=500, listw=listw_NY)
#summary(res)
#ecar1fw <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
# listw=listw_NY, weights=POP8, family="CAR", method="eigen")
#summary(ecar1fw)
#res <- MCMCsamp(ecar1fw, mcmc=5000, burnin=500, listw=listw_NY)
#summary(res)
esar0 <- errorsarlm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY)
summary(esar0)
res <- MCMCsamp(esar0, mcmc=5000, burnin=500, listw=listw_NY)
summary(res)
#esar0w <- errorsarlm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
# listw=listw_NY, weights=POP8)
#summary(esar0w)
#res <- MCMCsamp(esar0w, mcmc=5000, burnin=500, listw=listw_NY)
#summary(res)
esar1 <- errorsarlm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, etype="emixed")
summary(esar1)
res <- MCMCsamp(esar1, mcmc=5000, burnin=500, listw=listw_NY)
summary(res)
lsar0 <- lagsarlm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY)
summary(lsar0)
res <- MCMCsamp(lsar0, mcmc=5000, burnin=500, listw=listw_NY)
summary(res)
lsar1 <- lagsarlm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, type="mixed")
summary(lsar1)
res <- MCMCsamp(lsar1, mcmc=5000, burnin=500, listw=listw_NY)
summary(res)
ssar0 <- sacsarlm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY)
summary(ssar0)
res <- MCMCsamp(ssar0, mcmc=5000, burnin=500, listw=listw_NY)
summary(res)
ssar1 <- sacsarlm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, type="sacmixed")
summary(ssar1)
res <- MCMCsamp(ssar1, mcmc=5000, burnin=500, listw=listw_NY)
summary(res)
nydata <- st_read(system.file("shapes/NY8_bna_utm18.gpkg", package="spData")[1], quiet=TRUE)
lm0 <- lm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata)
summary(lm0)
lm0w <- lm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata, weights=POP8)
summary(lm0w)
suppressMessages(nyadjmat <- as.matrix(foreign::read.dbf(system.file(
 "misc/nyadjwts.dbf", package="spData")[1])[-1]))
suppressMessages(ID <- as.character(names(foreign::read.dbf(system.file(
 "misc/nyadjwts.dbf", package="spData")[1]))[-1]))
identical(substring(ID, 2, 10), substring(as.character(nydata$AREAKEY), 2, 10))
nyadjlw <- mat2listw(nyadjmat, as.character(nydata$AREAKEY))
listw_NY <- nb2listw(nyadjlw$neighbours, style="B")
esar0 <- errorsarlm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY)
summary(esar0)
system.time(esar1f <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
 data=nydata, listw=listw_NY, family="SAR", method="eigen"))
res <- summary(esar1f)
print(res)
sqrt(diag(res$resvar))
sqrt(diag(esar1f$fit$imat)*esar1f$fit$s2)
sqrt(diag(esar1f$fdHess))
system.time(esar1M <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
 data=nydata, listw=listw_NY, family="SAR", method="Matrix"))
summary(esar1M)
system.time(esar1M <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
 data=nydata, listw=listw_NY, family="SAR", method="Matrix",
 control=list(super=TRUE)))
summary(esar1M)
#esar1wf <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
# listw=listw_NY, weights=POP8, family="SAR", method="eigen")
#summary(esar1wf)
#system.time(esar1wM <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
# data=nydata, listw=listw_NY, weights=POP8, family="SAR", method="Matrix"))
#summary(esar1wM)
#esar1wlu <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
# listw=listw_NY, weights=POP8, family="SAR", method="LU")
#summary(esar1wlu)
#esar1wch <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
# listw=listw_NY, weights=POP8, family="SAR", method="Chebyshev")
#summary(esar1wch)
ecar1f <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, family="CAR", method="eigen")
summary(ecar1f)
system.time(ecar1M <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
 data=nydata, listw=listw_NY, family="CAR", method="Matrix"))
summary(ecar1M)
#ecar1wf <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
# listw=listw_NY, weights=nydata$POP8, family="CAR", method="eigen")
#summary(ecar1wf)
#system.time(ecar1wM <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
# data=nydata, listw=listw_NY, weights=POP8, family="CAR", method="Matrix"))
#summary(ecar1wM)
nc.sids <- st_read(system.file("shapes/sids.shp", package="spData")[1], quiet=TRUE)
ft.SID74 <- sqrt(1000)*(sqrt(nc.sids$SID74/nc.sids$BIR74) +
 sqrt((nc.sids$SID74+1)/nc.sids$BIR74))
lm_nc <- lm(ft.SID74 ~ 1)
sids.nhbr30 <- dnearneigh(cbind(nc.sids$east, nc.sids$north), 0, 30, row.names=row.names(nc.sids))
sids.nhbr30.dist <- nbdists(sids.nhbr30, cbind(nc.sids$east, nc.sids$north))
sids.nhbr <- listw2sn(nb2listw(sids.nhbr30, glist=sids.nhbr30.dist, style="B", zero.policy=TRUE))
dij <- sids.nhbr[,3]
n <- nc.sids$BIR74
el1 <- min(dij)/dij
el2 <- sqrt(n[sids.nhbr$to]/n[sids.nhbr$from])
sids.nhbr$weights <- el1*el2
sids.nhbr.listw <- sn2listw(sids.nhbr)
both <- factor(paste(nc.sids$L_id, nc.sids$M_id, sep=":"))
ft.NWBIR74 <- sqrt(1000)*(sqrt(nc.sids$NWBIR74/nc.sids$BIR74) +
 sqrt((nc.sids$NWBIR74+1)/nc.sids$BIR74))
mdata <- data.frame(both, ft.NWBIR74, ft.SID74, BIR74=nc.sids$BIR74)
outl <- which.max(rstandard(lm_nc))
as.character(nc.sids$NAME[outl])
mdata.4 <- mdata[-outl,]
W <- listw2mat(sids.nhbr.listw)
W.4 <- W[-outl, -outl]
sids.nhbr.listw.4 <- mat2listw(W.4)
esarI <- errorsarlm(ft.SID74 ~ 1, data=mdata, listw=sids.nhbr.listw,
 zero.policy=TRUE)
summary(esarI)
esarIa <- spautolm(ft.SID74 ~ 1, data=mdata, listw=sids.nhbr.listw,
 family="SAR")
summary(esarIa)
esarIV <- errorsarlm(ft.SID74 ~ ft.NWBIR74, data=mdata, listw=sids.nhbr.listw,
 zero.policy=TRUE)
summary(esarIV)
esarIVa <- spautolm(ft.SID74 ~ ft.NWBIR74, data=mdata, listw=sids.nhbr.listw,
 family="SAR")
summary(esarIVa)
#esarIaw <- spautolm(ft.SID74 ~ 1, data=mdata, listw=sids.nhbr.listw,
# weights=BIR74, family="SAR")
#summary(esarIaw)
#esarIIaw <- spautolm(ft.SID74 ~ both - 1, data=mdata, listw=sids.nhbr.listw,
# weights=BIR74, family="SAR")
#summary(esarIIaw)
#esarIVaw <- spautolm(ft.SID74 ~ ft.NWBIR74, data=mdata,
# listw=sids.nhbr.listw, weights=BIR74, family="SAR")
#summary(esarIVaw)
#ecarIaw <- spautolm(ft.SID74 ~ 1, data=mdata.4, listw=sids.nhbr.listw.4,
# weights=BIR74, family="CAR")
#summary(ecarIaw)
#ecarIIaw <- spautolm(ft.SID74 ~ both - 1, data=mdata.4,
# listw=sids.nhbr.listw.4, weights=BIR74, family="CAR")
#summary(ecarIIaw)
#ecarIVaw <- spautolm(ft.SID74 ~ ft.NWBIR74, data=mdata.4,
# listw=sids.nhbr.listw.4, weights=BIR74, family="CAR")
#summary(ecarIVaw)
#nc.sids$fitIV <- append(fitted.values(ecarIVaw), NA, outl-1)
#plot(nc.sids[,"fitIV"], nbreaks=12) # Cressie 1993, p. 565
data(oldcol)
COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="W"))
summary(COL.errW.eig)
COL.errW.sar <- spautolm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="W"))
summary(COL.errW.sar)
data(boston, package="spData")
gp1 <- spautolm(log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2)
 + I(RM^2) + AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT), 
 data=boston.c, nb2listw(boston.soi), family="SMA")
summary(gp1)
data(oldcol)
listw <- nb2listw(COL.nb, style="W")
ev <- eigenw(listw)
COL.sacW.eig <- sacsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, listw,
 control=list(pre_eig1=ev, pre_eig2=ev))
summary(COL.sacW.eig)
W <- as(listw, "CsparseMatrix")
trMatc <- trW(W, type="mult")
set.seed(1)
summary(impacts(COL.sacW.eig, tr=trMatc, R=2000), zstats=TRUE, short=TRUE)
library(coda)
set.seed(1)
COL.sacW.B0 <- spBreg_sac(CRIME ~ INC + HOVAL, data=COL.OLD, listw,
 Durbin=FALSE, control=list(ndraw=1500L, nomit=500L))
print(summary(COL.sacW.B0))
print(summary(impacts(COL.sacW.B0, tr=trMatc), zstats=TRUE, short=TRUE))
set.seed(1)
COL.sacW.B1 <- spBreg_sac(CRIME ~ INC + HOVAL, data=COL.OLD, listw,
 Durbin=TRUE, control=list(ndraw=1500L, nomit=500L))
print(summary(COL.sacW.B1))
print(summary(impacts(COL.sacW.B1, tr=trMatc), zstats=TRUE, short=TRUE))
COL.msacW.eig <- sacsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, listw,
 type="sacmixed", control=list(pre_eig1=ev, pre_eig2=ev))
summary(COL.msacW.eig)
set.seed(1)
summary(impacts(COL.msacW.eig, tr=trMatc, R=2000), zstats=TRUE, short=TRUE)
COL.msacW1.eig <- sacsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, listw,
 Durbin=TRUE, control=list(pre_eig1=ev, pre_eig2=ev))
summary(COL.msacW1.eig)
set.seed(1)
summary(impacts(COL.msacW1.eig, tr=trMatc, R=2000), zstats=TRUE, short=TRUE)
COL.msacW2.eig <- sacsarlm(CRIME ~ DISCBD + INC + HOVAL, data=COL.OLD, 
 listw, Durbin= ~ INC, control=list(pre_eig1=ev, pre_eig2=ev))
summary(COL.msacW2.eig)
summary(impacts(COL.msacW2.eig, tr=trMatc, R=2000), zstats=TRUE, short=TRUE)
library(coda)
set.seed(1)
data(oldcol)
lw <- nb2listw(COL.nb, style="W")
COL.err.Bayes <- spBreg_err(CRIME ~ INC + HOVAL, data=COL.OLD, listw=lw)
print(summary(COL.err.Bayes))
print(raftery.diag(COL.err.Bayes, r=0.01))
set.seed(1)
COL.err.Bayes <- spBreg_err(CRIME ~ INC + HOVAL, data=COL.OLD, listw=lw,
 control=list(prior=list(lambdaMH=TRUE)))
print(summary(COL.err.Bayes))
print(raftery.diag(COL.err.Bayes, r=0.01))
set.seed(1)
COL.err.Bayes <- spBreg_err(CRIME ~ INC + HOVAL, data=COL.OLD, listw=lw,
 Durbin=TRUE)
print(summary(COL.err.Bayes))
print(summary(impacts(COL.err.Bayes)))
print(raftery.diag(COL.err.Bayes, r=0.01))
set.seed(1)
COL.err.Bayes <- spBreg_err(CRIME ~ INC + HOVAL, data=COL.OLD, listw=lw,
 Durbin=TRUE, control=list(prior=list(lambdaMH=TRUE)))
print(summary(COL.err.Bayes))
print(summary(impacts(COL.err.Bayes)))
print(raftery.diag(COL.err.Bayes, r=0.01))
set.seed(1)
COL.err.Bayes <- spBreg_err(CRIME ~ INC + HOVAL, data=COL.OLD, listw=lw,
 Durbin=~INC)
print(summary(COL.err.Bayes))
print(summary(impacts(COL.err.Bayes)))
print(raftery.diag(COL.err.Bayes, r=0.01))
set.seed(1)
COL.err.Bayes <- spBreg_err(CRIME ~ INC + HOVAL, data=COL.OLD, listw=lw,
 Durbin=~INC, control=list(prior=list(lambdaMH=TRUE)))
print(summary(COL.err.Bayes))
print(summary(impacts(COL.err.Bayes)))
print(raftery.diag(COL.err.Bayes, r=0.01))
set.seed(1)
data(oldcol)
listw <- nb2listw(COL.nb, style="W")
COL.lag.Bayes <- spBreg_lag(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw=listw)
summary(COL.lag.Bayes)
summary(impacts(COL.lag.Bayes, tr=trMatc), short=TRUE, zstats=TRUE)
summary(impacts(COL.lag.Bayes, evalues=ev), short=TRUE, zstats=TRUE)
set.seed(1)
COL.D0.Bayes <- spBreg_lag(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw=listw, Durbin=TRUE)
summary(COL.D0.Bayes)
summary(impacts(COL.D0.Bayes, tr=trMatc), short=TRUE, zstats=TRUE)
set.seed(1)
COL.D1.Bayes <- spBreg_lag(CRIME ~ DISCBD + INC + HOVAL, data=COL.OLD,
 listw=listw, Durbin= ~ INC)
summary(COL.D1.Bayes)
summary(impacts(COL.D1.Bayes, tr=trMatc), short=TRUE, zstats=TRUE)
data(elect80, package="spData")
lw <- nb2listw(e80_queen, zero.policy=TRUE)
el_ml <- lagsarlm(log(pc_turnout) ~ log(pc_college) + log(pc_homeownership)
 + log(pc_income), data=elect80, listw=lw, zero.policy=TRUE, method="LU")
summary(el_ml)
set.seed(1)
el_B <- spBreg_lag(log(pc_turnout) ~ log(pc_college) + log(pc_homeownership)
 + log(pc_income), data=elect80, listw=lw, zero.policy=TRUE)
summary(el_B)
el_ml$timings
attr(el_B, "timings")
data(oldcol, package="spdep")
lw <- spdep::nb2listw(COL.nb)
COL.SLX <- lmSLX(CRIME ~ INC + HOVAL, data=COL.OLD, listw=lw)
pslx0 <- predict(COL.SLX)
pslx1 <- predict(COL.SLX, newdata=COL.OLD, listw=lw)
all.equal(pslx0, pslx1)
COL.OLD1 <- COL.OLD
COL.OLD1$INC <- COL.OLD1$INC + 1
pslx2 <- predict(COL.SLX, newdata=COL.OLD1, listw=lw)
sum(coef(COL.SLX)[c(2,4)])
mean(pslx2-pslx1)
COL.SLX <- lmSLX(CRIME ~ INC + HOVAL, data=COL.OLD, listw=lw)
summary(COL.SLX)
summary(impacts(COL.SLX))
COL.SLX <- lmSLX(CRIME ~ INC + HOVAL + I(HOVAL^2), data=COL.OLD, listw=lw, Durbin=TRUE)
summary(impacts(COL.SLX))
summary(COL.SLX)
COL.SLX <- lmSLX(CRIME ~ INC + HOVAL + I(HOVAL^2), data=COL.OLD, listw=lw, Durbin=~INC)
summary(impacts(COL.SLX))
summary(COL.SLX)
COL.SLX <- lmSLX(CRIME ~ INC, data=COL.OLD, listw=lw)
summary(COL.SLX)
summary(impacts(COL.SLX))
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
columbus <- st_read(system.file("shapes/columbus.shp", package="spData")[1], quiet=TRUE)
col.gal.nb <- read.gal(system.file("weights/columbus.gal", package="spData")[1])
# fit models for comparison
lm.mod <- lm(CRIME ~ HOVAL + INC, data=columbus)
lag <- lagsarlm(CRIME ~ HOVAL + INC, data=columbus, nb2listw(col.gal.nb))
mixed <- lagsarlm(CRIME ~ HOVAL + INC, data=columbus, nb2listw(col.gal.nb),
  Durbin=TRUE)
error <- errorsarlm(CRIME ~ HOVAL + INC, data=columbus, nb2listw(col.gal.nb))
# compare nested models
LR.sarlm(mixed, error)
#anova(lag, lm.mod)
#anova(lag, error, mixed)
AIC(lag, error, mixed)
columbus <- st_read(system.file("shapes/columbus.shp", package="spData")[1], quiet=TRUE)
col.gal.nb <- read.gal(system.file("weights/columbus.gal", package="spData")[1])
error.col <- errorsarlm(CRIME ~ HOVAL + INC, data=columbus,
 nb2listw(col.gal.nb))
bptest.sarlm(error.col)
bptest.sarlm(error.col, studentize=FALSE)
lm.target <- lm(error.col$tary ~ error.col$tarX - 1)
if (require(lmtest) && require(sandwich)) {
  print(coeftest(lm.target, vcov=vcovHC(lm.target, type="HC0"), df=Inf))
}
data(oldcol)
lw <- nb2listw(COL.nb, style="W")
ev <- eigenw(similar.listw(lw))
COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, quiet=FALSE, control=list(pre_eig=ev))
summary(COL.errW.eig)
COL.errW.eig_ev <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, control=list(pre_eig=ev))
all.equal(coefficients(COL.errW.eig), coefficients(COL.errW.eig_ev))
COL.errB.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="B"))
summary(COL.errB.eig)
W <- as(nb2listw(COL.nb), "CsparseMatrix")
trMatc <- trW(W, type="mult")
COL.errW.M <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="Matrix", quiet=FALSE, trs=trMatc)
summary(COL.errW.M)
COL.SDEM.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, etype="emixed", control=list(pre_eig=ev))
summary(COL.SDEM.eig)
COL.SDEM.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, Durbin=TRUE, control=list(pre_eig=ev))
summary(COL.SDEM.eig)
COL.SDEM.eig <- errorsarlm(CRIME ~ DISCBD + INC + HOVAL, data=COL.OLD,
 lw, Durbin=~INC, control=list(pre_eig=ev))
summary(COL.SDEM.eig)
summary(impacts(COL.SDEM.eig))
NA.COL.OLD <- COL.OLD
NA.COL.OLD$CRIME[20:25] <- NA
COL.err.NA <- errorsarlm(CRIME ~ INC + HOVAL, data=NA.COL.OLD,
 nb2listw(COL.nb), na.action=na.exclude)
COL.err.NA$na.action
COL.err.NA
resid(COL.err.NA)
lw <- nb2listw(COL.nb, style="W")
print(system.time(ev <- eigenw(similar.listw(lw))))
print(system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="eigen", control=list(pre_eig=ev))))
ocoef <- coefficients(COL.errW.eig)
print(system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="eigen", control=list(pre_eig=ev, LAPACK=FALSE))))
print(all.equal(ocoef, coefficients(COL.errW.eig)))
print(system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="eigen", control=list(pre_eig=ev, compiled_sse=TRUE))))
print(all.equal(ocoef, coefficients(COL.errW.eig)))
print(system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="Matrix_J", control=list(super=TRUE))))
print(all.equal(ocoef, coefficients(COL.errW.eig)))
print(system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="Matrix_J", control=list(super=FALSE))))
print(all.equal(ocoef, coefficients(COL.errW.eig)))
print(system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="Matrix_J", control=list(super=as.logical(NA)))))
print(all.equal(ocoef, coefficients(COL.errW.eig)))
print(system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="Matrix", control=list(super=TRUE))))
print(all.equal(ocoef, coefficients(COL.errW.eig)))
print(system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="Matrix", control=list(super=FALSE))))
print(all.equal(ocoef, coefficients(COL.errW.eig)))
print(system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="Matrix", control=list(super=as.logical(NA)))))
print(all.equal(ocoef, coefficients(COL.errW.eig)))
print(system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="spam", control=list(spamPivot="MMD"))))
print(all.equal(ocoef, coefficients(COL.errW.eig)))
print(system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="spam", control=list(spamPivot="RCM"))))
print(all.equal(ocoef, coefficients(COL.errW.eig)))
print(system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="spam_update", control=list(spamPivot="MMD"))))
print(all.equal(ocoef, coefficients(COL.errW.eig)))
print(system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 lw, method="spam_update", control=list(spamPivot="RCM"))))
print(all.equal(ocoef, coefficients(COL.errW.eig)))
columbus <- st_read(system.file("shapes/columbus.shp", package="spData")[1], quiet=TRUE)
col.gal.nb <- read.gal(system.file("weights/columbus.gal", package="spData")[1])
listw <- nb2listw(col.gal.nb)
ev <- eigenw(listw)
lobj <- lagsarlm(CRIME ~ INC + HOVAL, columbus, listw,
 control=list(pre_eig=ev))
summary(lobj)
mobj <- lagsarlm(CRIME ~ INC + HOVAL, columbus, listw, Durbin=TRUE,
 control=list(pre_eig=ev))
summary(mobj)
mobj1 <- lagsarlm(CRIME ~ INC + HOVAL, columbus, listw, Durbin= ~ INC,
 control=list(pre_eig=ev))
summary(mobj1)
W <- as(listw, "CsparseMatrix")
trMatc <- trW(W, type="mult")
trMC <- trW(W, type="MC")
set.seed(1)
impacts(lobj, listw=listw)
impacts(lobj, tr=trMatc)
impacts(lobj, tr=trMC)
impacts(lobj, evalues=ev)
library(coda)
lobjIQ5 <- impacts(lobj, tr=trMatc, R=200, Q=5)
summary(lobjIQ5, zstats=TRUE, short=TRUE)
summary(lobjIQ5, zstats=TRUE, short=TRUE, reportQ=TRUE)
impacts(mobj, listw=listw)
impacts(mobj, tr=trMatc)
impacts(mobj, tr=trMC)
impacts(mobj1, tr=trMatc)
impacts(mobj1, listw=listw)
cat(try(impacts(mobj, evalues=ev), silent=TRUE), "\n")
summary(impacts(mobj, tr=trMatc, R=200), short=TRUE, zstats=TRUE)
summary(impacts(mobj1, tr=trMatc, R=200), short=TRUE, zstats=TRUE)
#xobj <- lmSLX(CRIME ~ INC + HOVAL, columbus, listw)
#summary(impacts(xobj))
eobj <- errorsarlm(CRIME ~ INC + HOVAL, columbus, listw, etype="emixed")
summary(impacts(eobj), adjust_k=TRUE)
mobj1 <- lagsarlm(CRIME ~ INC + HOVAL, columbus, listw, type="mixed", 
method="Matrix", control=list(fdHess=TRUE))
summary(mobj1)
set.seed(1)
summary(impacts(mobj1, tr=trMatc, R=1000), zstats=TRUE, short=TRUE)
summary(impacts(mobj, tr=trMatc, R=1000), zstats=TRUE, short=TRUE)
mobj2 <- lagsarlm(CRIME ~ INC + HOVAL, columbus, listw, type="mixed", 
method="Matrix", control=list(fdHess=TRUE, optimHess=TRUE))
summary(impacts(mobj2, tr=trMatc, R=1000), zstats=TRUE, short=TRUE)
mobj3 <- lagsarlm(CRIME ~ INC + HOVAL, columbus, listw, type="mixed", 
method="spam", control=list(fdHess=TRUE))
summary(impacts(mobj3, tr=trMatc, R=1000), zstats=TRUE, short=TRUE)
data(boston, package="spData")
Wb <- as(nb2listw(boston.soi), "CsparseMatrix")
trMatb <- trW(Wb, type="mult")
gp2mMi <- lagsarlm(log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + 
I(RM^2) +  AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT), 
data=boston.c, nb2listw(boston.soi), type="mixed", method="Matrix", 
control=list(fdHess=TRUE), trs=trMatb)
summary(gp2mMi)
summary(impacts(gp2mMi, tr=trMatb, R=1000), zstats=TRUE, short=TRUE)
data(house, package="spData")
lw <- nb2listw(LO_nb)
form <- formula(log(price) ~ age + I(age^2) + I(age^3) + log(lotsize) +
   rooms + log(TLA) + beds + syear)
lobj <- lagsarlm(form, house, lw, method="Matrix",
 control=list(fdHess=TRUE), trs=trMat)
summary(lobj)
loobj <- impacts(lobj, tr=trMat, R=1000)
summary(loobj, zstats=TRUE, short=TRUE)
lobj1 <- stsls(form, house, lw)
loobj1 <- impacts(lobj1, tr=trMat, R=1000)
summary(loobj1, zstats=TRUE, short=TRUE)
mobj <- lagsarlm(form, house, lw, type="mixed",
 method="Matrix", control=list(fdHess=TRUE), trs=trMat)
summary(mobj)
moobj <- impacts(mobj, tr=trMat, R=1000)
summary(moobj, zstats=TRUE, short=TRUE)
data(oldcol)
listw <- nb2listw(COL.nb, style="W")
ev <- eigenw(listw)
W <- as(listw, "CsparseMatrix")
trMatc <- trW(W, type="mult")
COL.lag.eig <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, listw=listw,
 method="eigen", quiet=FALSE, control=list(pre_eig=ev, OrdVsign=1))
summary(COL.lag.eig, correlation=TRUE)
COL.lag.eig$fdHess
COL.lag.eig$resvar
# using the apparent sign in Ord (1975, equation B.1) 
COL.lag.eigb <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, listw=listw,
 method="eigen", control=list(pre_eig=ev, OrdVsign=-1))
summary(COL.lag.eigb)
COL.lag.eigb$fdHess
COL.lag.eigb$resvar
# force numerical Hessian
COL.lag.eig1 <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw=listw, method="Matrix", control=list(small=25))
summary(COL.lag.eig1)
COL.lag.eig1$fdHess
# force LeSage & Pace (2008, p. 57) approximation 
COL.lag.eig1a <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw=listw, method="Matrix", control=list(small=25), trs=trMatc)
summary(COL.lag.eig1a)
COL.lag.eig1a$fdHess
COL.lag.eig$resvar[2,2]
# using the apparent sign in Ord (1975, equation B.1) 
COL.lag.eigb$resvar[2,2]
# force numerical Hessian
COL.lag.eig1$fdHess[1,1]
# force LeSage & Pace (2008, p. 57) approximation 
COL.lag.eig1a$fdHess[2,2]
system.time(COL.lag.M <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb), method="Matrix", quiet=FALSE))
summary(COL.lag.M)
impacts(COL.lag.M, listw=nb2listw(COL.nb))
system.time(COL.lag.sp <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb), method="spam", quiet=FALSE))
summary(COL.lag.sp)
COL.lag.B <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="B"))
summary(COL.lag.B)
COL.mixed.B <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="B"), type="mixed", tol.solve=1e-9,
 control=list(pre_eig=ev))
summary(COL.mixed.B)
COL.mixed.W <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, type="mixed",
 control=list(pre_eig=ev))
summary(COL.mixed.W)
COL.mixed.D00 <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, Durbin=TRUE,
 control=list(pre_eig=ev))
summary(COL.mixed.D00)
COL.mixed.D01 <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, Durbin=FALSE,
 control=list(pre_eig=ev))
summary(COL.mixed.D01)
COL.mixed.D1 <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, Durbin= ~ INC + HOVAL,
 control=list(pre_eig=ev))
summary(COL.mixed.D1)
f <- CRIME ~ INC + HOVAL
COL.mixed.D2 <- lagsarlm(f, data=COL.OLD, listw,
 Durbin=as.formula(delete.response(terms(f))),
 control=list(pre_eig=ev))
summary(COL.mixed.D2)
COL.mixed.D1a <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, Durbin= ~ INC,
 control=list(pre_eig=ev))
summary(COL.mixed.D1a)
try(COL.mixed.D1 <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, Durbin= ~ inc + HOVAL,
 control=list(pre_eig=ev)))
try(COL.mixed.D1 <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 listw, Durbin= ~ DISCBD + HOVAL,
 control=list(pre_eig=ev)))
NA.COL.OLD <- COL.OLD
NA.COL.OLD$CRIME[20:25] <- NA
COL.lag.NA <- lagsarlm(CRIME ~ INC + HOVAL, data=NA.COL.OLD,
 nb2listw(COL.nb), na.action=na.exclude, 
 control=list(tol.opt=.Machine$double.eps^0.4))
COL.lag.NA$na.action
COL.lag.NA
resid(COL.lag.NA)
data(boston, package="spData")
gp2mM <- lagsarlm(log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + 
I(RM^2) +  AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT), 
data=boston.c, nb2listw(boston.soi), type="mixed", method="Matrix")
summary(gp2mM)
W <- as(nb2listw(boston.soi), "CsparseMatrix")
trMatb <- trW(W, type="mult")
gp2mMi <- lagsarlm(log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + 
I(RM^2) +  AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT), 
data=boston.c, nb2listw(boston.soi), type="mixed", method="Matrix", 
trs=trMatb)
summary(gp2mMi)
data(oldcol)
lw <- nb2listw(COL.nb)
COL.lag.eig <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, lw)

COL.mix.eig <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, lw,
  type="mixed")
print(p1 <- predict(COL.mix.eig))
#print(p2 <- predict(COL.mix.eig, newdata=COL.OLD, listw=lw, pred.type = "TS",
# legacy.mixed = TRUE))
AIC(COL.mix.eig)
sqrt(deviance(COL.mix.eig)/length(COL.nb))
sqrt(sum((COL.OLD$CRIME - as.vector(p1))^2)/length(COL.nb))
#sqrt(sum((COL.OLD$CRIME - as.vector(p2))^2)/length(COL.nb))

COL.err.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, lw)
AIC(COL.err.eig)
sqrt(deviance(COL.err.eig)/length(COL.nb))
sqrt(sum((COL.OLD$CRIME - as.vector(predict(COL.err.eig)))^2)/length(COL.nb))
#sqrt(sum((COL.OLD$CRIME - as.vector(predict(COL.err.eig, newdata=COL.OLD,
#  listw=lw, pred.type = "TS")))^2)/length(COL.nb))

COL.SDerr.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, lw,
 etype="emixed")
AIC(COL.SDerr.eig)
sqrt(deviance(COL.SDerr.eig)/length(COL.nb))
sqrt(sum((COL.OLD$CRIME - as.vector(predict(COL.SDerr.eig)))^2)/length(COL.nb))
#sqrt(sum((COL.OLD$CRIME - as.vector(predict(COL.SDerr.eig, newdata=COL.OLD,
#  listw=lw, pred.type = "TS")))^2)/length(COL.nb))

AIC(COL.lag.eig)
sqrt(deviance(COL.lag.eig)/length(COL.nb))
sqrt(sum((COL.OLD$CRIME - as.vector(predict(COL.lag.eig)))^2)/length(COL.nb))
#sqrt(sum((COL.OLD$CRIME - as.vector(predict(COL.lag.eig, newdata=COL.OLD,
#  listw=lw, pred.type = "TS")))^2)/length(COL.nb))

#p3 <- predict(COL.mix.eig, newdata=COL.OLD, listw=lw, pred.type = "TS",
# legacy=FALSE, legacy.mixed = TRUE)
#all.equal(p2, p3, check.attributes=FALSE)
#p4 <- predict(COL.mix.eig, newdata=COL.OLD, listw=lw, pred.type = "TS",
# legacy=FALSE, power=TRUE, legacy.mixed = TRUE)
#all.equal(p2, p4, check.attributes=FALSE)
#p5 <- predict(COL.mix.eig, newdata=COL.OLD, listw=lw, pred.type = "TS",
# legacy=TRUE, power=TRUE, legacy.mixed = TRUE)
#all.equal(p2, p5, check.attributes=FALSE)
data(oldcol)
listw <- nb2listw(COL.nb, style="W")
ev <- eigenw(listw)
COL.sacW.eig <- sacsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, listw,
 control=list(pre_eig1=ev, pre_eig2=ev))
summary(COL.sacW.eig)
W <- as(listw, "CsparseMatrix")
trMatc <- trW(W, type="mult")
set.seed(1)
summary(impacts(COL.sacW.eig, tr=trMatc, R=2000), zstats=TRUE, short=TRUE)
COL.msacW.eig <- sacsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, listw,
 type="sacmixed", control=list(pre_eig1=ev, pre_eig2=ev))
summary(COL.msacW.eig)
set.seed(1)
summary(impacts(COL.msacW.eig, tr=trMatc, R=2000), zstats=TRUE, short=TRUE)
COL.msacW1.eig <- sacsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, listw,
 Durbin=TRUE, control=list(pre_eig1=ev, pre_eig2=ev))
summary(COL.msacW1.eig)
set.seed(1)
summary(impacts(COL.msacW1.eig, tr=trMatc, R=2000), zstats=TRUE, short=TRUE)
COL.msacW2.eig <- sacsarlm(CRIME ~ DISCBD + INC + HOVAL, data=COL.OLD, 
 listw, Durbin= ~ INC, control=list(pre_eig1=ev, pre_eig2=ev))
summary(COL.msacW2.eig)
summary(impacts(COL.msacW2.eig, tr=trMatc, R=2000), zstats=TRUE, short=TRUE)
data(oldcol)
COL.mix.eig <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb), type="mixed", method="eigen")
summary(COL.mix.eig, correlation=TRUE, Nagelkerke=TRUE)
COL.mix.M <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb), type="mixed", method="Matrix")
summary(COL.mix.M, correlation=TRUE, Nagelkerke=TRUE)
COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
  nb2listw(COL.nb, style="W"), method="eigen")
summary(COL.errW.eig, correlation=TRUE, Nagelkerke=TRUE, Hausman=TRUE)

} #dontrun
}


\keyword{spatial}

