\name{spdep-defunct}
\alias{spdep-defunct}
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
\alias{impacts.gmsar}
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
\alias{trW}


\title{Defunct Functions in Package \pkg{spdep}}

\description{
  These functions are defunct from release 1.2-1. The
  functions have been moved to the \pkg{spatialreg} package.
}

\usage{
lextrB(lw, zero.policy = TRUE, control = list())
lextrW(lw, zero.policy=TRUE, control=list())
lextrS(lw, zero.policy=TRUE, control=list())
#l_max(lw, zero.policy=TRUE, control=list())
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
trW(W=NULL, m = 30, p = 16, type = "mult", listw=NULL, momentsSymmetry=TRUE)
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
  \item{momentsSymmetry}{default TRUE; assert Smirnov/Anselin symmetry assumption}
}

\details{
 Model-fitting functions and functions supporting model fitting have been moved to the \pkg{spatialreg} package.
}


\seealso{
  \code{\link[base]{Defunct}}
}


\keyword{spatial}

