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


\title{Deprecated Functions in Package \pkg{spdep}}

\description{
  These functions are provided for compatibility with older versions of
  \pkg{spdep} only, and may be defunct as soon as the next release. The
  functions have been moved to the \pkg{spreg} package.
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

}
\section{Control arguments: lextrB/W/S, l_max}{
\describe{
  \item{trace}{report values in while loops, default NULL assuming FALSE; logical}
  \item{tol}{tolerance for breaking while loops, default \code{.Machine$double.eps^(1/2)}; numeric}
  \item{maxiter}{maximum number of iterations in while loops, default \code{6 * (length(lw$neighbours) - 2}; integer}
  \item{useC}{use C code, default TRUE, logical (not in \code{l_max})}
}
}

\details{
 Model-fitting functions and functions supporting model fitting are being moved to the \pkg{spreg} package.
The \code{griffith_sone} function function may be used, following Ord and Gasim (for references see Griffith and Sone (1995)), to calculate analytical eigenvalues for binary rook or queen contiguous neighbours where the data are arranged as a regular P times Q grid. The \code{subgraph_eigenw} function may be used when there are multiple graph components, of which the largest may be handled as a dense matrix. Here the eigenvalues are computed for each subgraph in turn, and catenated to reconstruct the complete set. The functions may be used to provide pre-computed eigenvalues for spatial regression functions.
}


\seealso{
  \code{\link[base]{Deprecated}}
}

\examples{
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
\dontrun{
ab.listw <- nb2listw(boston.soi, style="S")
er <- range(eigenw(similar.listw(ab.listw)))
er
res_1 <- lextrS(ab.listw)
c(res_1)
}
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
\donttest{
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
}
\dontrun{
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
}
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

}

\keyword{spatial}

