\name{impacts}
\Rdversion{1.1}
\alias{impacts}
\alias{impacts.sarlm}
\alias{impacts.stsls}
\alias{impacts.gmsar}
\alias{impacts.SLX}
\alias{impacts.MCMC_sar_g}
\alias{impacts.lagmess}
\alias{plot.lagImpact}
\alias{print.lagImpact}
\alias{summary.lagImpact}
\alias{print.summary.lagImpact}
\alias{print.WXImpact}
\alias{summary.WXImpact}
\alias{print.summary.WXImpact}
\alias{HPDinterval.lagImpact}
\alias{intImpacts}
\title{Impacts in spatial lag models}
\description{The calculation of impacts for spatial lag and spatial Durbin models is needed in order to interpret the regression coefficients correctly, because of the spillovers between the terms in these data generation processes (unlike the spatial error model). Methods for \dQuote{SLX} and Bayesian fitted models are also provided, the former do not need MC simulations, while the latter pass through MCMC draws.}
\usage{
\method{impacts}{sarlm}(obj, \dots, tr, R = NULL, listw = NULL, useHESS = NULL,
 tol = 1e-06, empirical = FALSE, Q=NULL)
\method{impacts}{stsls}(obj, \dots, tr, R = NULL, listw = NULL,
 tol = 1e-06, empirical = FALSE, Q=NULL)
\method{impacts}{gmsar}(obj, \dots, n = NULL, tr = NULL, R = NULL, listw = NULL,
 tol = 1e-06, empirical = FALSE, Q=NULL)
\method{impacts}{lagmess}(obj, ..., R=NULL, listw=NULL, tol=1e-6, empirical=FALSE)
\method{impacts}{SLX}(obj, ...)
\method{impacts}{MCMC_sar_g}(obj, ..., tr=NULL, listw=NULL, Q=NULL)
\method{plot}{lagImpact}(x, ..., choice="direct", trace=FALSE, density=TRUE)
\method{print}{lagImpact}(x, ..., reportQ=NULL)
\method{summary}{lagImpact}(object, ..., zstats=FALSE, short=FALSE, reportQ=NULL)
\method{print}{WXImpact}(x, ...)
\method{summary}{WXImpact}(object, ..., adjust_k=FALSE)
\method{HPDinterval}{lagImpact}(obj, prob = 0.95, ..., choice="direct")
intImpacts(rho, beta, P, n, mu, Sigma, irho, drop2beta, bnames, interval,
 type, tr, R, listw, tol, empirical, Q, icept, iicept, p, mess=FALSE,
 samples=NULL)
}
\arguments{
  \item{obj}{A spatial regression object created by \code{lagsarlm}, \code{lagmess} or by \code{lmSLX}; in \code{HPDinterval.lagImpact}, a lagImpact object}
  \item{\dots}{Arguments passed through to methods in the \pkg{coda} package}
  \item{tr}{A vector of traces of powers of the spatial weights matrix created using \code{trW}, for approximate impact measures; if not given, \code{listw} must be given for exact measures (for small to moderate spatial weights matrices); the traces must be for the same spatial weights as were used in fitting the spatial regression, and must be row-standardised}
  \item{listw}{If \code{tr} is not given, a spatial weights object as created by \code{nb2listw}; they must be the same spatial weights as were used in fitting the spatial regression, but do not have to be row-standardised}
  \item{n}{defaults to \code{length(obj$residuals)}; in the method for \code{gmsar} objects it may be used in panel settings to compute the impacts for cross-sectional weights only, suggested by Angela Parenti}
  \item{R}{If given, simulations are used to compute distributions for the impact measures, returned as \code{mcmc} objects; the objects are used for convenience but are not output by an MCMC process}
  \item{useHESS}{Use the Hessian approximation (if available) even if the asymptotic coefficient covariance matrix is available; used for comparing methods}
  \item{tol}{Argument passed to \code{mvrnorm}: tolerance (relative to largest variance) for numerical lack of positive-definiteness in the coefficient covariance matrix}
  \item{empirical}{Argument passed to \code{mvrnorm} (default FALSE): if true, the coefficients and their covariance matrix specify the empirical not population mean and covariance matrix}
  \item{Q}{default NULL, else an integer number of cumulative power series impacts to calculate if \code{tr} is given}
  \item{reportQ}{default NULL; if TRUE and \code{Q} given as an argument to \code{impacts}, report impact components}
  \item{x, object}{lagImpact objects created by \code{impacts} methods}
  \item{zstats}{default FALSE, if TRUE, also return z-values and p-values for the impacts based on the simulations}
  \item{short}{default FALSE, if TRUE passed to the print summary method to omit printing of the mcmc summaries}
  \item{choice}{One of three impacts: direct, indirect, or total}
  \item{trace}{Argument passed to \code{plot.mcmc}: plot trace plots}
  \item{density}{Argument passed to \code{plot.mcmc}: plot density plots}
  \item{prob}{Argument passed to \code{HPDinterval.mcmc}: a numeric scalar in the interval (0,1) giving the target probability content of the intervals}
  \item{adjust_k}{adjust ML SDEM standard errors by dividing by (n-k) rather than n}
  \item{rho, beta, P, mu, Sigma, irho, drop2beta, bnames,
    interval, type, icept, iicept, p, mess, samples}{internal arguments shared inside impacts methods}
}

\details{If called without \code{R} being set, the method returns the direct, indirect and total impacts for the variables in the model, for the variables themselves in tha spatial lag model case, for the variables and their spatial lags in the spatial Durbin (mixed) model case. The spatial lag impact measures are computed using eq. 2.46 (LeSage and Pace, 2009, p. 38), either using the exact dense matrix (when \code{listw} is given), or traces of powers of the weights matrix (when \code{tr} is given). When the traces are created by powering sparse matrices, the exact and the trace methods should give very similar results, unless the number of powers used is very small, or the spatial coefficient is close to its bounds.

If \code{R} is given, simulations will be used to create distributions for the impact measures, provided that the fitted model object contains a coefficient covariance matrix. The simulations are made using \code{\link[MASS]{mvrnorm}} with the coefficients and their covariance matrix from the fitted model.

The simulations are stored as \code{mcmc} objects as defined in the \pkg{coda} package; the objects are used for convenience but are not output by an MCMC process. The simulated values of the coefficients are checked to see that the spatial coefficient remains within its valid interval --- draws outside the interval are discarded.

When \code{Q} and \code{tr} are given, addition impact component results are provided for each step in the traces of powers of the weights matrix up to and including the \code{Q}'th power. This increases computing time because the output object is substantially increased in size in proportion to the size of \code{Q}.

The method for \code{gmsar} objects is only for those of \code{type} \code{SARAR} output by \code{gstsls}, and assume that the spatial error coefficient is fixed, and thus omitted from the coefficients and covariance matrix used for simulation.
}

\value{
An object of class lagImpact.

If no simulation is carried out, the object returned is a list with:
\item{direct}{numeric vector}
\item{indirect}{numeric vector}
\item{total}{numeric vector}
and a matching \code{Qres} list attribute if \code{Q} was given.

If simulation is carried out, the object returned is a list with:

\item{res}{a list with three components as for the non-simulation case, with a matching \code{Qres} list attribute if \code{Q} was given}
\item{sres}{a list with three \code{mcmc} matrices, for the direct, indirect and total impacts with a matching \code{Qmcmc} list attribute if \code{Q} was given}
}

\references{LeSage J and RK Pace (2009) \emph{Introduction to Spatial Econometrics}. CRC Press, Boca Raton, pp. 33--42, 114--115; LeSage J and MM Fischer (2008) Spatial growth regressions: model specification, estimation and interpretation. \emph{Spatial Economic Analysis} 3 (3), pp. 275--304.

Roger Bivand, Gianfranco Piras (2015). Comparing Implementations of Estimation Methods for Spatial Econometrics. \emph{Journal of Statistical Software}, 63(18), 1-36. \url{https://www.jstatsoft.org/v63/i18/}.
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}


\seealso{
\code{\link{trW}}, \code{\link{lagsarlm}}, \code{\link{nb2listw}}, \code{\link[MASS]{mvrnorm}}, \code{\link[coda]{plot.mcmc}}, \code{\link[coda]{summary.mcmc}}, \code{\link[coda]{HPDinterval}}
}
\examples{
if (require(rgdal, quietly=TRUE)) {
example(columbus, package="spData")
listw <- nb2listw(col.gal.nb)
lobj <- lagsarlm(CRIME ~ INC + HOVAL, columbus, listw)
summary(lobj)
mobj <- lagsarlm(CRIME ~ INC + HOVAL, columbus, listw, type="mixed")
summary(mobj)
W <- as(listw, "CsparseMatrix")
trMatc <- trW(W, type="mult")
trMC <- trW(W, type="MC")
set.seed(1)
impacts(lobj, listw=listw)
impacts(lobj, tr=trMatc)
impacts(lobj, tr=trMC)
lobj1 <- stsls(CRIME ~ INC + HOVAL, columbus, listw)
loobj1 <- impacts(lobj1, tr=trMatc, R=200)
summary(loobj1, zstats=TRUE, short=TRUE)
library(coda)
HPDinterval(loobj1)
lobj1r <- stsls(CRIME ~ INC + HOVAL, columbus, listw, robust=TRUE)
loobj1r <- impacts(lobj1r, tr=trMatc, R=200)
summary(loobj1r, zstats=TRUE, short=TRUE)
lobjIQ5 <- impacts(lobj, tr=trMatc, R=200, Q=5)
summary(lobjIQ5, zstats=TRUE, short=TRUE)
summary(lobjIQ5, zstats=TRUE, short=TRUE, reportQ=TRUE)
impacts(mobj, listw=listw)
impacts(mobj, tr=trMatc)
impacts(mobj, tr=trMC)
summary(impacts(mobj, tr=trMatc, R=200), zstats=TRUE)
xobj <- lmSLX(CRIME ~ INC + HOVAL, columbus, listw)
summary(impacts(xobj))
eobj <- errorsarlm(CRIME ~ INC + HOVAL, columbus, listw, etype="emixed")
summary(impacts(eobj), adjust_k=TRUE)
\dontrun{
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
}
\dontrun{
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
}
}
}
\keyword{spatial}

