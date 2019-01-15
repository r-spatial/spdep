\name{aple}
\alias{aple}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Approximate profile-likelihood estimator (APLE)}
\description{
  The Approximate profile-likelihood estimator (APLE) of the simultaneous autoregressive model's spatial dependence parameter was introduced in Li et al. (2007). It employs a correction term using the eigenvalues of the spatial weights matrix, and consequently should not be used for large numbers of observations. It also requires that the variable has a mean of zero, and it is assumed that it has been detrended. The spatial weights object is assumed to be row-standardised, that is using default \code{style="W"} in \code{nb2listw}.
}
\usage{
aple(x, listw, override_similarity_check=FALSE, useTrace=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a zero-mean detrended continuous variable}
  \item{listw}{a \code{listw} object from for example \code{nb2listw}}
  \item{override\_similarity\_check}{default FALSE, if TRUE - typically for row-standardised weights with asymmetric underlying general weights - similarity is not checked}
  \item{useTrace}{default TRUE, use trace of sparse matrix \code{W \%*\% W} (Li et al. (2010)), if FALSE, use crossproduct of eigenvalues of \code{W} as in Li et al. (2007)}
}
\details{
 This implementation has been checked with Hongfei Li's own implementation using her data; her help was very valuable.
}
\value{
  A scalar APLE value.
}
\references{Li, H, Calder, C. A. and Cressie N. A. C. (2007) Beyond Moran's I: testing for spatial dependence based on the spatial autoregressive model. Geographical Analysis 39, 357-375; Li, H, Calder, C. A. and Cressie N. A. C. (2012) One-step estimation of spatial dependence parameters: Properties and extensions of the APLE statistic, Journal of Multivariate Analysis 105, 68-84.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{nb2listw}}, \code{\link{aple.mc}}, \code{\link{aple.plot}}
}
\examples{
run <- FALSE
if (require(sf, quietly=TRUE)) run <- TRUE
if (run) {
wheat <- st_read(system.file("shapes/wheat.shp", package="spData")[1], quiet=TRUE)
nbr1 <- poly2nb(as(wheat, "Spatial"), queen=FALSE)
nbrl <- nblag(nbr1, 2)
nbr12 <- nblag_cumul(nbrl)
wheat_g <- wheat
st_geometry(wheat_g) <- NULL
cms0 <- with(wheat_g, tapply(yield, c, median))
cms1 <- c(model.matrix(~ factor(c) -1, data=wheat) \%*\% cms0)
wheat$yield_detrend <- wheat$yield - cms1
wheat_g$yield_detrend <- wheat_g$yield - cms1
isTRUE(all.equal(c(with(wheat_g,
 tapply(yield_detrend, c, median))), rep(0.0, 25),
 check.attributes=FALSE))
}
if (run) {
moran.test(wheat$yield_detrend, nb2listw(nbr12, style="W"))
}
if (run) {
aple(as.vector(scale(wheat$yield_detrend, scale=FALSE)), nb2listw(nbr12, style="W"))
}
\dontrun{
if (run) {
errorsarlm(yield_detrend ~ 1, wheat, nb2listw(nbr12, style="W"))
}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}

