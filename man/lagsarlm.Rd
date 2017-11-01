% Copyright 2002-16 by Roger S. Bivand
\name{lagsarlm}
\alias{lagsarlm}
\alias{spBreg_lag}
\title{Spatial simultaneous autoregressive lag model estimation}
\description{
  The \code{lagsarlm} function provides Maximum likelihood estimation of spatial simultaneous autoregressive lag and spatial Durbin (mixed) models of the form:

\deqn{y = \rho W y + X \beta + \varepsilon}{y = rho W y + X beta + e}

where \eqn{\rho}{rho} is found by \code{optimize()} first, and \eqn{\beta}{beta} and other parameters by generalized least squares subsequently (one-dimensional search using optim performs badly on some platforms). In the spatial Durbin (mixed) model, the spatially lagged independent variables are added to X. Note that interpretation of the fitted coefficients should use impact measures, because of the feedback loops induced by the data generation process for this model. With one of the sparse matrix methods, larger numbers of observations can be handled, but the \code{interval=} argument may need be set when the weights are not row-standardised.

The \code{spBreg_lag} function is an early-release version of the Matlab Spatial Econometrics Toolbox function \code{sar_g.m}, using drawing by inversion, and not accommodating heteroskedastic disturbances.
}
\usage{
lagsarlm(formula, data = list(), listw, 
	na.action, type="lag", method="eigen", quiet=NULL, 
	zero.policy=NULL, interval=NULL, tol.solve=1.0e-10, trs=NULL, 
	control=list())
spBreg_lag(formula, data = list(), listw, na.action, type="lag",
    zero.policy=NULL, control=list())
}
\arguments{
  \item{formula}{a symbolic description of the model to be fit. The details 
of model specification are given for \code{lm()}}
  \item{data}{an optional data frame containing the variables in the model. 
By default the variables are taken from the environment which the function 
is called.}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{na.action}{a function (default \code{options("na.action")}), can also be \code{na.omit} or \code{na.exclude} with consequences for residuals and fitted values - in these cases the weights list will be subsetted to remove NAs in the data. It may be necessary to set zero.policy to TRUE because this subsetting may create no-neighbour observations. Note that only weights lists created without using the glist argument to \code{nb2listw} may be subsetted.}
  \item{type}{default "lag", may be set to "mixed"; when "mixed", the lagged intercept is dropped for spatial weights style "W", that is row-standardised weights, but otherwise included; \dQuote{Durbin} may be used instead of \dQuote{mixed}}
  \item{method}{"eigen" (default) - the Jacobian is computed as the product 
of (1 - rho*eigenvalue) using \code{eigenw}, and "spam" or "Matrix_J" for strictly symmetric weights lists of styles "B" and "C", or made symmetric by similarity (Ord, 1975, Appendix C) if possible for styles "W" and "S", using code from the spam or Matrix packages to calculate the determinant; \dQuote{Matrix} and \dQuote{spam_update} provide updating Cholesky decomposition methods; "LU" provides an alternative sparse matrix decomposition approach. In addition, there are "Chebyshev" and Monte Carlo "MC" approximate log-determinant methods; the Smirnov/Anselin (2009) trace approximation is available as "moments". Three methods: "SE_classic", "SE_whichMin", and "SE_interp" are provided experimentally, the first to attempt to emulate the behaviour of Spatial Econometrics toolbox ML fitting functions. All use grids of log determinant values, and the latter two attempt to ameliorate some features of "SE_classic".
}
  \item{quiet}{default NULL, use !verbose global option value; if FALSE, reports function values during optimization.}
  \item{zero.policy}{default NULL, use global option value; if TRUE assign zero to the lagged value of zones without 
neighbours, if FALSE (default) assign NA - causing \code{lagsarlm()} to terminate with an error}
  \item{interval}{default is NULL, search interval for autoregressive parameter}
  \item{tol.solve}{the tolerance for detecting linear dependencies in the columns of matrices to be inverted - passed to \code{solve()} (default=1.0e-10). This may be used if necessary to extract coefficient standard errors (for instance lowering to 1e-12), but errors in \code{solve()} may constitute indications of poorly scaled variables: if the variables have scales differing much from the autoregressive coefficient, the values in this matrix may be very different in scale, and inverting such a matrix is analytically possible by definition, but numerically unstable; rescaling the RHS variables alleviates this better than setting tol.solve to a very small value}
  \item{trs}{default NULL, if given, a vector of powered spatial weights matrix traces output by \code{trW}; when given, insert the asymptotic analytical values into the numerical Hessian instead of the approximated values; may be used to get around some problems raised when the numerical Hessian is poorly conditioned, generating NaNs in subsequent operations; the use of trs is recommended}
 \item{control}{list of extra control arguments - see section below}
}
\details{
  The asymptotic standard error of \eqn{\rho}{rho} is only computed when
method=eigen, because the full matrix operations involved would be costly
for large n typically associated with the choice of method="spam" or "Matrix". The same applies to the coefficient covariance matrix. Taken as the
asymptotic matrix from the literature, it is typically badly scaled, and with the elements involving \eqn{\rho}{rho} being very small,
while other parts of the matrix can be very large (often many orders
of magnitude in difference). It often happens that the \code{tol.solve}
argument needs to be set to a smaller value than the default, or the RHS variables can be centred or reduced in range.

Versions of the package from 0.4-38 include numerical Hessian values where asymptotic standard errors are not available. This change has been introduced to permit the simulation of distributions for impact measures. The warnings made above with regard to variable scaling also apply in this case.

Note that the fitted() function for the output object assumes that the response 
variable may be reconstructed as the sum of the trend, the signal, and the
noise (residuals). Since the values of the response variable are known,
their spatial lags are used to calculate signal components (Cressie 1993, p. 564). This differs from other software, including GeoDa, which does not use knowledge of the response 
variable in making predictions for the fitting data.
}

\section{Control arguments}{
\describe{
  \item{tol.opt:}{the desired accuracy of the optimization - passed to \code{optimize()} (default=square root of double precision machine tolerance, a larger root may be used needed, see help(boston) for an example)}
  \item{fdHess:}{default NULL, then set to (method != "eigen") internally; use \code{fdHess} to compute an approximate Hessian using finite differences when using sparse matrix methods; used to make a coefficient covariance matrix when the number of observations is large; may be turned off to save resources if need be}
  \item{optimHess:}{default FALSE, use \code{fdHess} from \pkg{nlme}, if TRUE, use \code{optim} to calculate Hessian at optimum}
  \item{optimHessMethod:}{default \dQuote{optimHess}, may be \dQuote{nlm} or one of the \code{optim} methods}
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

\section{Extra Bayesian control arguments}{
\describe{
  \item{ldet_method}{default \dQuote{SE_classic}; equivalent to the \code{method} argument in \code{lagsarlm}}
  \item{interval}{default \code{c(-1, 1)}; used unmodified or set internally by \code{jacobianSetup}}
  \item{ndraw}{default \code{2500L}; integer total number of draws}
  \item{nomit}{default \code{500L}; integer total number of omitted burn-in draws}
  \item{thin}{default \code{1L}; integer thinning proportion}
  \item{verbose}{default \code{FALSE}; inverse of \code{quiet} argument in \code{lagsarlm}}
  \item{detval}{default \code{NULL}; not yet in use, precomputed matrix of log determinants}
  \item{prior}{a list with the following components:
    \describe{
      \item{Tbeta}{default \code{NULL}; values of the betas variance-covariance matrix, set to \code{diag(k)*1e+12} if \code{NULL}}
      \item{c_beta}{default \code{NULL}; values of the betas set to 0 if \code{NULL}}
      \item{rho}{default \code{0.5}; value of the autoregressive coefficient}
      \item{sige}{default \code{1}; value of the residual variance}
      \item{nu}{default \code{0}; informative Gamma(nu,d0) prior on sige}
      \item{d0}{default \code{0}; informative Gamma(nu,d0) prior on sige}
      \item{a1}{default \code{1.01}; parameter for beta(a1,a2) prior on rho}
      \item{a2}{default \code{1.01}; parameter for beta(a1,a2) prior on rho}
    }}
}}

\value{
  A list object of class \code{sarlm}
  \item{type}{"lag" or "mixed"}
  \item{rho}{simultaneous autoregressive lag coefficient}
  \item{coefficients}{GLS coefficient estimates}
  \item{rest.se}{asymptotic standard errors if ase=TRUE, otherwise approximate numeriacal Hessian-based values}
  \item{LL}{log likelihood value at computed optimum}
  \item{s2}{GLS residual variance}
  \item{SSE}{sum of squared GLS errors}
  \item{parameters}{number of parameters estimated}
%  \item{lm.model}{the \code{lm} object returned when estimating for \eqn{\rho=0}{rho=0}}
  \item{logLik_lm.model}{Log likelihood of the linear model for \eqn{\rho=0}{rho=0}}
  \item{AIC_lm.model}{AIC of the linear model for \eqn{\rho=0}{rho=0}}
  \item{method}{the method used to calculate the Jacobian}
  \item{call}{the call used to create this object}
  \item{residuals}{GLS residuals}
  \item{tarX}{model matrix of the GLS model}
  \item{tary}{response of the GLS model}
  \item{y}{response of the linear model for \eqn{\rho=0}{rho=0}}
  \item{X}{model matrix of the linear model for \eqn{\rho=0}{rho=0}}
%  \item{lm.target}{the \code{lm} object returned for the GLS fit}
  \item{opt}{object returned from numerical optimisation}
  \item{fitted.values}{Difference between residuals and response variable}
  \item{se.fit}{Not used yet}
%  \item{formula}{model formula}
  \item{ase}{TRUE if method=eigen}
  \item{rho.se}{if ase=TRUE, the asymptotic standard error of \eqn{\rho}{rho}, otherwise approximate numeriacal Hessian-based value}
  \item{LMtest}{if ase=TRUE, the Lagrange Multiplier test for the absence
of spatial autocorrelation in the lag model residuals}
  \item{resvar}{the asymptotic coefficient covariance matrix for (s2, rho, B)}
  \item{zero.policy}{zero.policy for this model}
  \item{aliased}{the aliased explanatory variables (if any)}
  \item{listw_style}{the style of the spatial weights used}
  \item{interval}{the line search interval used to find \eqn{\rho}{rho}}
  \item{fdHess}{the numerical Hessian-based coefficient covariance matrix for (rho, B) if computed}
  \item{optimHess}{if TRUE and fdHess returned, \code{optim} used to calculate Hessian at optimum}
  \item{insert}{if TRUE and fdHess returned, the asymptotic analytical values are inserted into the numerical Hessian instead of the approximated values, and its size increased to include the first row/column for sigma2}
  \item{LLNullLlm}{Log-likelihood of the null linear model}
  \item{timings}{processing timings}
  \item{f_calls}{number of calls to the log likelihood function during optimization}
  \item{hf_calls}{number of calls to the log likelihood function during numerical Hessian computation}
  \item{intern_classic}{a data frame of detval matrix row choices used by the SE toolbox classic method}
  \item{na.action}{(possibly) named vector of excluded or omitted observations if non-default na.action argument used}

The internal sar.lag.mixed.* functions return the value of the log likelihood function at \eqn{\rho}{rho}.

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
\author{Roger Bivand \email{Roger.Bivand@nhh.no}, with thanks to Andrew 
Bernat for contributions to the asymptotic standard error code.}

\seealso{\code{\link{lm}}, \code{\link{errorsarlm}}, 
\code{\link{summary.sarlm}}, \code{\link{eigenw}},  
\code{\link{predict.sarlm}}, \code{\link{impacts.sarlm}},
\code{\link{residuals.sarlm}}, \code{\link{do_ldet}}
}

\examples{
data(oldcol)
COL.lag.eig <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="W"), method="eigen", quiet=FALSE)
summary(COL.lag.eig, correlation=TRUE)
COL.lag.eig$fdHess
COL.lag.eig$resvar
W <- as(nb2listw(COL.nb), "CsparseMatrix")
trMatc <- trW(W, type="mult")
COL.lag.eig1 <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="W"), control=list(fdHess=TRUE), trs=trMatc)
COL.lag.eig1$fdHess
system.time(COL.lag.M <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb), method="Matrix", quiet=FALSE))
summary(COL.lag.M)
impacts(COL.lag.M, listw=nb2listw(COL.nb))
\dontrun{
system.time(COL.lag.sp <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb), method="spam", quiet=FALSE))
summary(COL.lag.sp)
}
COL.lag.B <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="B"))
summary(COL.lag.B, correlation=TRUE)
COL.mixed.B <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="B"), type="mixed", tol.solve=1e-9)
summary(COL.mixed.B, correlation=TRUE)
COL.mixed.W <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="W"), type="mixed")
summary(COL.mixed.W, correlation=TRUE)
NA.COL.OLD <- COL.OLD
NA.COL.OLD$CRIME[20:25] <- NA
COL.lag.NA <- lagsarlm(CRIME ~ INC + HOVAL, data=NA.COL.OLD,
 nb2listw(COL.nb), na.action=na.exclude, 
 control=list(tol.opt=.Machine$double.eps^0.4))
COL.lag.NA$na.action
COL.lag.NA
resid(COL.lag.NA)
\dontrun{
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
}
summary(COL.lag.eig)
COL.lag.Bayes <- spBreg_lag(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="W"))
summary(COL.lag.Bayes)
set.seed(1)
summary(impacts(COL.lag.Bayes, tr=trMatc), short=TRUE, zstats=TRUE)
\dontrun{
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
}
}
\keyword{spatial}
