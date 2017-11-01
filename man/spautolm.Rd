\name{spautolm}
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
%\alias{CAR}
%\alias{SAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Spatial conditional and simultaneous autoregression model estimation}
\description{
 Function taking family and weights arguments for spatial autoregression model estimation by Maximum Likelihood, using dense matrix methods, not suited to large data sets with thousands of observations. With one of the sparse matrix methods, larger numbers of observations can be handled, but the \code{interval=} argument should be set. The implementation is GLS using the single spatial coefficient value, here termed lambda, found by line search using \code{optimize} to maximise the log likelihood.
}
\usage{
spautolm(formula, data = list(), listw, weights,
 na.action, family = "SAR", method="eigen", verbose = NULL, trs=NULL,
 interval=NULL, zero.policy = NULL, tol.solve=.Machine$double.eps,
 llprof=NULL, control=list())
\method{summary}{spautolm}(object, correlation = FALSE, adj.se=FALSE,
 Nagelkerke=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a symbolic description of the model to be fit. The details 
of model specification are given for \code{lm()}}
  \item{data}{an optional data frame containing the variables in the model. 
By default the variables are taken from the environment which the function 
is called.}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{weights}{an optional vector of weights to be used in the fitting process}
  \item{na.action}{a function (default \code{options("na.action")}), can also be \code{na.omit} or \code{na.exclude} with consequences for residuals and fitted values - in these cases the weights list will be subsetted to remove NAs in the data. Note that only weights lists created without using the glist argument to \code{nb2listw} may be subsetted.}
  \item{family}{character string: either \code{"SAR"} or \code{"CAR"} for simultaneous or conditional autoregressions; \code{"SMA"} for spatial moving average added thanks to Jielai Ma - \code{"SMA"} is only implemented for method=\code{"eigen"} because it necessarily involves dense matrices}
  \item{method}{character string: default \code{"eigen"} for use of dense matrices, \code{"Matrix_J"} for sparse matrices (restricted to spatial weights symmetric or similar to symmetric) using methods in the Matrix package; \dQuote{Matrix} provides updating Cholesky decomposition methods. Values of method may also include "LU", which provides an alternative sparse matrix decomposition approach, and the "Chebyshev" and Monte Carlo "MC" approximate log-determinant methods.}
  \item{verbose}{default NULL, use global option value; if TRUE, reports function values during optimization.}
  \item{trs,}{default NULL, if given, a vector of powered spatial weights matrix traces output by \code{trW}; when given, used in some Jacobian methods}
  \item{interval}{search interval for autoregressive parameter when not using method="eigen"; default is c(-1,0.999), \code{optimize} will reset NA/NaN to a bound and gives a warning when the interval is poorly set; method="Matrix" will attempt to search for an appropriate interval, if find\_interval=TRUE (fails on some platforms)}
  \item{zero.policy}{default NULL, use global option value; Include list of no-neighbour observations in output if TRUE --- otherwise zero.policy is handled within the listw argument}
  \item{tol.solve}{the tolerance for detecting linear dependencies in the columns of matrices to be inverted - passed to \code{solve()} (default=double precision machine tolerance). Errors in \code{solve()} may constitute indications of poorly scaled variables: if the variables have scales differing much from the autoregressive coefficient, the values in this matrix may be very different in scale, and inverting such a matrix is analytically possible by definition, but numerically unstable; rescaling the RHS variables alleviates this better than setting tol.solve to a very small value}
  \item{llprof}{default NULL, can either be an integer, to divide the feasible range into llprof points, or a sequence of spatial coefficient values, at which to evaluate the likelihood function}
 \item{control}{list of extra control arguments - see section below}
  \item{object}{\code{spautolm} object from \code{spautolm}}
  \item{correlation}{logical; if 'TRUE', the correlation matrix of the estimated parameters is returned and printed (default=FALSE)}
  \item{adj.se}{if TRUE, adjust the coefficient standard errors for the number of fitted coefficients}
  \item{Nagelkerke}{if TRUE, the Nagelkerke pseudo R-squared is reported}
  \item{\dots}{further arguments passed to or from other methods}
}


\details{
  This implementation is based on \code{\link[MASS]{lm.gls}} and \code{\link{errorsarlm}}. In particular, the function does not (yet) prevent asymmetric spatial weights being used with "CAR" family models. It appears that both numerical issues (convergence in particular) and uncertainties about the exact spatial weights matrix used make it difficult to reproduce Cressie and Chan's 1989 results, also given in Cressie 1993.

Note that the fitted() function for the output object assumes that the response 
variable may be reconstructed as the sum of the trend, the signal, and the
noise (residuals). Since the values of the response variable are known,
their spatial lags are used to calculate signal components (Cressie 1993, p. 564). This differs from other software, including GeoDa, which does not use knowledge of the response 
variable in making predictions for the fitting data.

}


\section{Control arguments}{
\describe{
  \item{tol.opt:}{the desired accuracy of the optimization - passed to \code{optimize()} (default=\code{.Machine$double.eps^(2/3)})}
  \item{fdHess:}{default NULL, then set to (method != "eigen") internally; use \code{fdHess} to compute an approximate Hessian using finite differences when using sparse matrix methods; used to make a coefficient covariance matrix when the number of observations is large; may be turned off to save resources if need be}
  \item{optimHess:}{default FALSE, use \code{fdHess} from \pkg{nlme}, if TRUE, use \code{optim} to calculate Hessian at optimum}
  \item{optimHessMethod:}{default \dQuote{optimHess}, may be \dQuote{nlm} or one of the \code{optim} methods}
  \item{Imult:}{default 2; used for preparing the Cholesky decompositions for updating in the Jacobian function}
  \item{super:}{if NULL (default), set to FALSE to use a simplicial decomposition for the sparse Cholesky decomposition and method \dQuote{Matrix_J}, set to  \code{as.logical(NA)} for method \dQuote{Matrix}, if TRUE, use a supernodal decomposition}
  \item{cheb_q:}{default 5; highest power of the approximating polynomial for the Chebyshev approximation}
  \item{MC_p:}{default 16; number of random variates}
  \item{MC_m:}{default 30; number of products of random variates matrix and spatial weights matrix}
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
  A list object of class \code{spautolm}:
  \item{fit}{a list, with items:
  \describe{
    \item{coefficients}{ML coefficient estimates}
    \item{SSE}{ML sum of squared errors}
    \item{s2}{ML residual variance}
    \item{imat}{ML coefficient covariance matrix (before multiplying by s2)}
    \item{signal\_trend}{non-spatial component of fitted.values}
    \item{signal\_stochastic}{spatial component of fitted.values}
    \item{fitted.values}{sum of non-spatial and spatial components of fitted.values}
    \item{residuals}{difference between observed and fitted values}
  }}
  \item{lambda}{ML autoregressive coefficient}
  \item{LL}{log likelihood for fitted model}
  \item{LL0}{log likelihood for model with lambda=0}
  \item{call}{the call used to create this object}
  \item{parameters}{number of parameters estimated}
  \item{aliased}{if not NULL, details of aliased variables}
  \item{method}{Jacobian method chosen}
  \item{family}{family chosen}
  \item{zero.policy}{zero.policy used}
  \item{weights}{case weights used}
  \item{interval}{the line search interval used}
  \item{timings}{processing timings}
  \item{na.action}{(possibly) named vector of excluded or omitted observations if non-default na.action argument used}
  \item{llprof}{if not NULL, a list with components lambda and ll of equal length}
  \item{lambda.se}{Numerical Hessian-based standard error of lambda}
  \item{fdHess}{Numerical Hessian-based variance-covariance matrix}
  \item{X}{covariates used in model fitting}
  \item{Y}{response used in model fitting}
  \item{weights}{weights used in model fitting}
}
\references{Cliff, A. D., Ord, J. K. 1981 \emph{Spatial processes}, Pion;
Ord, J. K. 1975 Estimation methods for models of spatial interaction,
\emph{Journal of the American Statistical Association}, 70, 120-126; Waller, L. A., Gotway, C. A. 2004 \emph{Applied spatial statistics for public health}, Wiley, Hoboken, NJ, 325-380; Cressie, N. A. C. 1993 \emph{Statistics for spatial data}, Wiley, New York, 548-568; Ripley, B. D. 1981 \emph{Spatial statistics}, Wiley, New York, 88-95; LeSage J and RK Pace (2009) Introduction to Spatial Econometrics. CRC Press, Boca Raton.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}
\note{The standard errors given in Waller and Gotway (2004) are adjusted for the numbers of parameters estimated, and may be reproduced by using the additional argument \code{adj.se=TRUE} in the \code{summary} method. In addition, the function returns fitted values and residuals as given by Cressie (1993) p. 564.}

\seealso{\code{\link[stats]{optimize}}, \code{\link{errorsarlm}}, \code{\link{do_ldet}}}
\examples{
\dontrun{
if (require(foreign, quietly=TRUE)) {
example(NY_data, package="spData")
lm0 <- lm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata)
summary(lm0)
lm0w <- lm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata, weights=POP8)
summary(lm0w)
esar0 <- errorsarlm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY)
summary(esar0)
system.time(esar1f <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
 data=nydata, listw=listw_NY, family="SAR", method="eigen", verbose=TRUE))
res <- summary(esar1f)
print(res)
sqrt(diag(res$resvar))
sqrt(diag(esar1f$fit$imat)*esar1f$fit$s2)
sqrt(diag(esar1f$fdHess))
system.time(esar1M <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
 data=nydata, listw=listw_NY, family="SAR", method="Matrix", verbose=TRUE))
summary(esar1M)
system.time(esar1M <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
 data=nydata, listw=listw_NY, family="SAR", method="Matrix", verbose=TRUE,
 control=list(super=TRUE)))
summary(esar1M)
esar1wf <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, weights=POP8, family="SAR", method="eigen")
summary(esar1wf)
system.time(esar1wM <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
 data=nydata, listw=listw_NY, weights=POP8, family="SAR", method="Matrix"))
summary(esar1wM)
esar1wlu <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, weights=POP8, family="SAR", method="LU")
summary(esar1wlu)
esar1wch <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, weights=POP8, family="SAR", method="Chebyshev")
summary(esar1wch)
ecar1f <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, family="CAR", method="eigen")
summary(ecar1f)
system.time(ecar1M <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
 data=nydata, listw=listw_NY, family="CAR", method="Matrix"))
summary(ecar1M)
ecar1wf <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, weights=nydata$POP8, family="CAR", method="eigen")
summary(ecar1wf)
system.time(ecar1wM <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
 data=nydata, listw=listw_NY, weights=POP8, family="CAR", method="Matrix"))
summary(ecar1wM)
}
if (require(rgdal, quietly=TRUE)) {
example(nc.sids, package="spData")
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
as.character(nc.sids$names[outl])
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
esarIaw <- spautolm(ft.SID74 ~ 1, data=mdata, listw=sids.nhbr.listw,
 weights=BIR74, family="SAR")
summary(esarIaw)
esarIIaw <- spautolm(ft.SID74 ~ both - 1, data=mdata, listw=sids.nhbr.listw,
 weights=BIR74, family="SAR")
summary(esarIIaw)
esarIVaw <- spautolm(ft.SID74 ~ ft.NWBIR74, data=mdata,
 listw=sids.nhbr.listw, weights=BIR74, family="SAR")
summary(esarIVaw)
ecarIaw <- spautolm(ft.SID74 ~ 1, data=mdata.4, listw=sids.nhbr.listw.4,
 weights=BIR74, family="CAR")
summary(ecarIaw)
ecarIIaw <- spautolm(ft.SID74 ~ both - 1, data=mdata.4,
 listw=sids.nhbr.listw.4, weights=BIR74, family="CAR")
summary(ecarIIaw)
ecarIVaw <- spautolm(ft.SID74 ~ ft.NWBIR74, data=mdata.4,
 listw=sids.nhbr.listw.4, weights=BIR74, family="CAR")
summary(ecarIVaw)
nc.sids$fitIV <- append(fitted.values(ecarIVaw), NA, outl-1)
spplot(nc.sids, c("fitIV"), cuts=12) # Cressie 1993, p. 565
}
\dontrun{
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
}
}
}
\keyword{spatial}

