\name{MCMCsamp}
\alias{MCMCsamp}
\alias{MCMCsamp.spautolm}
\alias{MCMCsamp.sarlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{MCMC sample from fitted spatial regression}
\description{The \code{MCMCsamp} method uses \code{\link[LearnBayes]{rwmetrop}}, a random walk Metropolis algorithm, from \pkg{LearnBayes} to make MCMC samples from fitted maximum likelihood spatial regression models.}
\usage{
MCMCsamp(object, mcmc = 1L, verbose = NULL, ...)
\method{MCMCsamp}{spautolm}(object, mcmc = 1L, verbose = NULL, ...,
 burnin = 0L, scale=1, listw, control = list())
\method{MCMCsamp}{sarlm}(object, mcmc = 1L, verbose = NULL, ...,
    burnin=0L, scale=1, listw, listw2=NULL, control=list())}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A spatial regression model object fitted by maximum likelihood with \code{\link{spautolm}}}
  \item{mcmc}{The number of MCMC iterations after burnin}
  \item{verbose}{default NULL, use global option value; if TRUE, reports progress}
  \item{\dots}{Arguments passed through}
  \item{burnin}{The number of burn-in iterations for the sampler}
  \item{scale}{a positive scale parameter}
  \item{listw, listw2}{\code{listw} objects created for example by \code{nb2listw}; should be the same object(s) used for fitting the model}
  \item{control}{list of extra control arguments - see \code{\link{spautolm}}}
}
\value{An object of class \dQuote{mcmc} suited to \pkg{coda}, with attributes: \dQuote{accept} acceptance rate; \dQuote{type} input ML fitted model type \dQuote{SAR}, \dQuote{CAR}, \dQuote{SMA}, \dQuote{lag}, \dQuote{mixed}, \dQuote{error}, \dQuote{sac}, \dQuote{sacmixed}; \dQuote{timings} run times}
\note{If the acceptance rate is below 0.05, a warning will be issued; consider increasing mcmc.}
\references{Jim Albert (2007) Bayesian Computation with R, Springer, New York, pp. 104-105.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link[LearnBayes]{rwmetrop}}, \code{\link{spautolm}}, \code{\link{lagsarlm}}, \code{\link{errorsarlm}}, \code{\link{sacsarlm}}}
\examples{
if (require(foreign, quietly=TRUE)) {
example(NY_data, package="spData")
\dontrun{
esar1f <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, family="SAR", method="eigen")
summary(esar1f)
res <- MCMCsamp(esar1f, mcmc=5000, burnin=500, listw=listw_NY)
summary(res)
esar1fw <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, weights=POP8, family="SAR", method="eigen")
summary(esar1fw)
res <- MCMCsamp(esar1fw, mcmc=5000, burnin=500, listw=listw_NY)
summary(res)
ecar1f <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, family="CAR", method="eigen")
summary(ecar1f)
res <- MCMCsamp(ecar1f, mcmc=5000, burnin=500, listw=listw_NY)
summary(res)
esar1fw <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, weights=POP8, family="SAR", method="eigen")
summary(esar1fw)
res <- MCMCsamp(esar1fw, mcmc=5000, burnin=500, listw=listw_NY)
summary(res)
ecar1fw <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, weights=POP8, family="CAR", method="eigen")
summary(ecar1fw)
res <- MCMCsamp(ecar1fw, mcmc=5000, burnin=500, listw=listw_NY)
summary(res)
}
esar0 <- errorsarlm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY)
summary(esar0)
res <- MCMCsamp(esar0, mcmc=5000, burnin=500, listw=listw_NY)
summary(res)
\dontrun{
esar0w <- errorsarlm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, weights=POP8)
summary(esar0)
res <- MCMCsamp(esar0w, mcmc=5000, burnin=500, listw=listw_NY)
summary(res)
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
}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}
