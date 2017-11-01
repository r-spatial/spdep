% Copyright 2010 by Roger S. Bivand
\name{sacsarlm}
\alias{sacsarlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Spatial simultaneous autoregressive SAC model estimation}
\description{
  Maximum likelihood estimation of spatial simultaneous autoregressive
\dQuote{SAC/SARAR} models of the form:

\deqn{y = \rho W1 y + X \beta + u, u = \lambda W2 u + \varepsilon}{y = rho W1 y + X beta + u, u = lambda W2 u + e}

where \eqn{\rho}{rho} and \eqn{\lambda}{lambda} are found by \code{nlminb} or \code{optim()} first, and \eqn{\beta}{beta} and other parameters by generalized least squares subsequently}
\usage{
sacsarlm(formula, data = list(), listw, listw2 = NULL, na.action, type="sac",
 method = "eigen", quiet = NULL, zero.policy = NULL, tol.solve = 1e-10,
 llprof=NULL, interval1=NULL, interval2=NULL, trs1=NULL, trs2=NULL,
 control = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a symbolic description of the model to be fit. The details 
of model specification are given for \code{lm()}}
  \item{data}{an optional data frame containing the variables in the model. 
By default the variables are taken from the environment which the function 
is called}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{listw2}{a \code{listw} object created for example by \code{nb2listw}, if not given, set to the same spatial weights as the \code{listw} argument}
  \item{na.action}{a function (default \code{options("na.action")}), can also be \code{na.omit} or \code{na.exclude} with consequences for residuals and fitted values - in these cases the weights list will be subsetted to remove NAs in the data. It may be necessary to set zero.policy to TRUE because this subsetting may create no-neighbour observations. Note that only weights lists created without using the glist argument to \code{nb2listw} may be subsetted.}
  \item{type}{default "sac", may be set to "sacmixed" for the Manski model to include the spatially lagged independent variables added to X using \code{listw}; when "sacmixed", the lagged intercept is dropped for spatial weights style "W", that is row-standardised weights, but otherwise included}
  \item{method}{"eigen" (default) - the Jacobian is computed as the product 
of (1 - rho*eigenvalue) using \code{eigenw}, and "spam" or "Matrix" for strictly symmetric weights lists of styles "B" and "C", or made symmetric by similarity (Ord, 1975, Appendix C) if possible for styles "W" and "S", using code from the spam or Matrix packages to calculate the determinant; "LU" provides an alternative sparse matrix decomposition approach. In addition, there are "Chebyshev" and Monte Carlo "MC" approximate log-determinant methods.}
  \item{quiet}{default NULL, use !verbose global option value; if FALSE, reports function values during optimization.}
  \item{zero.policy}{default NULL, use global option value; if TRUE assign zero to the lagged value of zones without 
neighbours, if FALSE (default) assign NA - causing \code{sacsarlm()} to terminate with an error}
  \item{tol.solve}{the tolerance for detecting linear dependencies in the columns of matrices to be inverted - passed to \code{solve()} (default=1.0e-10). This may be used if necessary to extract coefficient standard errors (for instance lowering to 1e-12), but errors in \code{solve()} may constitute indications of poorly scaled variables: if the variables have scales differing much from the autoregressive coefficient, the values in this matrix may be very different in scale, and inverting such a matrix is analytically possible by definition, but numerically unstable; rescaling the RHS variables alleviates this better than setting tol.solve to a very small value}
  \item{llprof}{default NULL, can either be an integer, to divide the feasible ranges into a grid of points, or a two-column matrix of spatial coefficient values, at which to evaluate the likelihood function}
  \item{trs1, trs2}{default NULL, if given, vectors for each weights object of powered spatial weights matrix traces output by \code{trW}; when given, used in some Jacobian methods}
  \item{interval1, interval2}{default is NULL, search intervals for each weights object for autoregressive parameters}
  \item{control}{list of extra control arguments - see section below}
}

\details{Because numerical optimisation is used to find the values of lambda and rho, care needs to be shown. It has been found that the surface of the 2D likelihood function often forms a \dQuote{banana trench} from (low rho, high lambda) through (high rho, high lambda) to (high rho, low lambda) values. In addition, sometimes the banana has optima towards both ends, one local, the other global, and conseqently the choice of the starting point for the final optimization becomes crucial. The default approach is not to use just (0, 0) as a starting point, nor the (rho, lambda) values from \code{gstsls}, which lie in a central part of the \dQuote{trench}, but either four values at (low rho, high lambda), (0, 0), (high rho, high lambda), and (high rho, low lambda), and to use the best of these start points for the final optimization. Optionally, nine points can be used spanning the whole (lower, upper) space.}

\section{Control arguments}{
\describe{
  \item{fdHess:}{default NULL, then set to (method != "eigen") internally; use \code{fdHess} to compute an approximate Hessian using finite differences when using sparse matrix methods with \code{fdHess} from \pkg{nlme}; used to make a coefficient covariance matrix when the number of observations is large; may be turned off to save resources if need be}
  \item{LAPACK:}{default FALSE; logical value passed to \code{qr} in the SSE log likelihood function}
  \item{Imult:}{default 2; used for preparing the Cholesky decompositions for updating in the Jacobian function}
  \item{cheb_q:}{default 5; highest power of the approximating polynomial for the Chebyshev approximation}
  \item{MC_p:}{default 16; number of random variates}
  \item{MC_m:}{default 30; number of products of random variates matrix and spatial weights matrix}
  \item{super:}{default FALSE using a simplicial decomposition for the sparse Cholesky decomposition, if TRUE, use a supernodal decomposition}
  \item{opt_method:}{default \dQuote{nlminb}, may be set to \dQuote{L-BFGS-B} to use box-constrained optimisation in \code{optim}}
  \item{opt_control:}{default \code{list()}, a control list to pass to \code{nlminb} or \code{optim}}
  \item{pars:}{default \code{NULL}, for which five trial starting values spanning the lower/upper range are tried and the best selected, starting values of \eqn{\rho}{rho} and \eqn{\lambda}{lambda}}
  \item{npars}{default integer \code{4L}, four trial points; if not default value, nine trial points}
  \item{pre_eig1, pre_eig2}{default NULL; may be used to pass pre-computed vectors of eigenvalues}
}}


\value{
  A list object of class \code{sarlm}
  \item{type}{\dQuote{sac}}
  \item{rho}{lag simultaneous autoregressive lag coefficient}
  \item{lambda}{error simultaneous autoregressive error coefficient}
  \item{coefficients}{GLS coefficient estimates}
  \item{rest.se}{asymptotic standard errors if ase=TRUE, otherwise approximate numeriacal Hessian-based values}
  \item{ase}{TRUE if method=eigen}
  \item{LL}{log likelihood value at computed optimum}
  \item{s2}{GLS residual variance}
  \item{SSE}{sum of squared GLS errors}
  \item{parameters}{number of parameters estimated}
  \item{logLik_lm.model}{Log likelihood of the non-spatial linear model}
  \item{AIC_lm.model}{AIC of the non-spatial linear model}
%  \item{lm.model}{the \code{lm} object returned when estimating for \eqn{\rho=0}{rho=0}}
  \item{method}{the method used to calculate the Jacobian}
  \item{call}{the call used to create this object}
  \item{residuals}{GLS residuals}
%  \item{lm.target}{the \code{lm} object returned for the GLS fit}
  \item{tarX}{model matrix of the GLS model}
  \item{tary}{response of the GLS model}
  \item{y}{response of the linear model for \eqn{\rho=0}{rho=0}}
  \item{X}{model matrix of the linear model for \eqn{\rho=0}{rho=0}}
  \item{opt}{object returned from numerical optimisation}
  \item{pars}{starting parameter values for final optimization, either given or found by trial point evaluation}
  \item{mxs}{if default input pars, optimal objective function values at trial points}
  \item{fitted.values}{Difference between residuals and response variable}
  \item{se.fit}{Not used yet}
%  \item{formula}{model formula}
  \item{rho.se}{if ase=TRUE, the asymptotic standard error of \eqn{\rho}{rho}, otherwise approximate numeriacal Hessian-based value}
  \item{lambda.se}{if ase=TRUE, the asymptotic standard error of \eqn{\lambda}{lambda}}
  \item{resvar}{the asymptotic coefficient covariance matrix for (s2, rho, lambda, B)}
  \item{zero.policy}{zero.policy for this model}
  \item{aliased}{the aliased explanatory variables (if any)}
  \item{LLNullLlm}{Log-likelihood of the null linear model}
  \item{fdHess}{the numerical Hessian-based coefficient covariance matrix for (rho, lambda, B) if computed}
  \item{resvar}{asymptotic coefficient covariance matrix}
  \item{optimHess}{FALSE}
  \item{timings}{processing timings}
  \item{na.action}{(possibly) named vector of excluded or omitted observations if non-default na.action argument used}
}
\references{Anselin, L. 1988 \emph{Spatial econometrics: methods and models.}
(Dordrecht: Kluwer); LeSage J and RK Pace (2009) \emph{Introduction to Spatial Econometrics}. CRC Press, Boca Raton.

Roger Bivand, Gianfranco Piras (2015). Comparing Implementations of Estimation Methods for Spatial Econometrics. \emph{Journal of Statistical Software}, 63(18), 1-36. \url{https://www.jstatsoft.org/v63/i18/}.

Bivand, R. S., Hauke, J., and Kossowski, T. (2013). Computing the Jacobian in Gaussian spatial autoregressive models: An illustrated comparison of available methods. \emph{Geographical Analysis}, 45(2), 150-179.
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{lm}}, \code{\link{lagsarlm}}, \code{\link{errorsarlm}}, 
\code{\link{summary.sarlm}}, \code{\link{eigenw}}, \code{\link{impacts.sarlm}}}
\examples{
data(oldcol)
COL.sacW.eig <- sacsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, 
 nb2listw(COL.nb, style="W"))
summary(COL.sacW.eig, correlation=TRUE)
W <- as(nb2listw(COL.nb, style="W"), "CsparseMatrix")
trMatc <- trW(W, type="mult")
summary(impacts(COL.sacW.eig, tr=trMatc, R=2000), zstats=TRUE, short=TRUE)
COL.msacW.eig <- sacsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, 
 nb2listw(COL.nb, style="W"), type="sacmixed")
summary(COL.msacW.eig, correlation=TRUE)
summary(impacts(COL.msacW.eig, tr=trMatc, R=2000), zstats=TRUE, short=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}
