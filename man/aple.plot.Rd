\name{aple.plot}
\alias{aple.plot}
\alias{localAple}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Approximate profile-likelihood estimator (APLE) scatterplot}
\description{
  A scatterplot decomposition of the approximate profile-likelihood estimator, and a local APLE based on the list of vectors returned by the scatterplot function.
}
\usage{
aple.plot(x, listw, override_similarity_check=FALSE, useTrace=TRUE, do.plot=TRUE, ...)
localAple(x, listw, override_similarity_check=FALSE, useTrace=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a zero-mean detrended continuous variable}
  \item{listw}{a \code{listw} object from for example \code{nb2listw}}
  \item{override\_similarity\_check}{default FALSE, if TRUE - typically for row-standardised weights with asymmetric underlying general weights - similarity is not checked}
  \item{useTrace}{default TRUE, use trace of sparse matrix \code{W \%*\% W} (Li et al. (2010)), if FALSE, use crossproduct of eigenvalues of \code{W} as in Li et al. (2007)}
  \item{do.plot}{default TRUE: should a scatterplot be drawn}
  \item{\dots}{other arguments to be passed to \code{plot}}
}
\details{
  The function solves a secondary eigenproblem of size n internally, so constructing the values for the scatterplot is quite compute and memory intensive, and is not suitable for very large n.
}
\value{
\code{aple.plot} returns list with components:
  \item{X}{A vector as described in Li et al. (2007), p. 366.}
  \item{Y}{A vector as described in Li et al. (2007), p. 367.}

\code{localAple} returns a vector of local APLE values.
}
\references{Li, H, Calder, C. A. and Cressie N. A. C. (2007) Beyond Moran's I: testing for spatial dependence based on the spatial autoregressive model. Geographical Analysis 39, pp. 357-375; Li, H, Calder, C. A. and Cressie N. A. C. (2012) One-step estimation of spatial dependence parameters: Properties and extensions of the APLE statistic, Journal of Multivariate Analysis 105, 68-84.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{aple}}}
\examples{
\dontrun{
if (require(rgdal, quietly=TRUE)) {
example(aple)
plt_out <- aple.plot(as.vector(scale(wheat$yield_detrend, scale=FALSE)),
 nb2listw(nbr12, style="W"), cex=0.6)
crossprod(plt_out$Y, plt_out$X)/crossprod(plt_out$X)
lm_obj <- lm(Y ~ X, plt_out)
abline(lm_obj)
abline(v=0, h=0, lty=2)
zz <- summary(influence.measures(lm_obj))
infl <- as.integer(rownames(zz))
points(plt_out$X[infl], plt_out$Y[infl], pch=3, cex=0.6, col="red")
wheat$localAple <- localAple(as.vector(scale(wheat$yield_detrend, scale=FALSE)),
 nb2listw(nbr12, style="W"))
mean(wheat$localAple)
hist(wheat$localAple)
spl <- list("sp.text", coordinates(wheat)[infl,], rep("*", length(infl)))
spplot(wheat, "localAple", sp.layout=spl)
}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}

