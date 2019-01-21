\name{aple.mc}
\alias{aple.mc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Approximate profile-likelihood estimator (APLE) permutation test}
\description{
  A permutation bootstrap test for the approximate profile-likelihood estimator (APLE).
}
\usage{
aple.mc(x, listw, nsim, override_similarity_check=FALSE, useTrace=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a zero-mean detrended continuous variable}
  \item{listw}{a \code{listw} object from for example \code{nb2listw}}
  \item{nsim}{number of simulations}
  \item{override\_similarity\_check}{default FALSE, if TRUE - typically for row-standardised weights with asymmetric underlying general weights - similarity is not checked}
  \item{useTrace}{default TRUE, use trace of sparse matrix \code{W \%*\% W} (Li et al. (2010)), if FALSE, use crossproduct of eigenvalues of \code{W} as in Li et al. (2007)}
}

\value{
  A \code{boot} object as returned by the \code{boot} function.
}
\references{Li, H, Calder, C. A. and Cressie N. A. C. (2007) Beyond Moran's I: testing for spatial dependence based on the spatial autoregressive model. Geographical Analysis 39, 357-375; Li, H, Calder, C. A. and Cressie N. A. C. (2012) One-step estimation of spatial dependence parameters: Properties and extensions of the APLE statistic, Journal of Multivariate Analysis 105, 68-84.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{aple}}, \code{\link[boot]{boot}}}
\examples{
\dontrun{
wheat <- st_read(system.file("shapes/wheat.shp", package="spData")[1], quiet=TRUE)
nbr1 <- poly2nb(wheat, queen=FALSE)
nbrl <- nblag(nbr1, 2)
nbr12 <- nblag_cumul(nbrl)
wheat_g <- wheat
st_geometry(wheat_g) <- NULL
cms0 <- with(wheat_g, tapply(yield, c, median))
cms1 <- c(model.matrix(~ factor(c) -1, data=wheat) \%*\% cms0)
wheat$yield_detrend <- wheat$yield - cms1
oldRNG <- RNGkind()
RNGkind("L'Ecuyer-CMRG")
set.seed(1L)
boot_out_ser <- aple.mc(as.vector(scale(wheat$yield_detrend, scale=FALSE)),
 nb2listw(nbr12, style="W"), nsim=500)
plot(boot_out_ser)
boot_out_ser
library(parallel)
oldCores <- set.coresOption(NULL)
nc <- detectCores(logical=FALSE)
# set nc to 1L here
if (nc > 1L) nc <- 1L
invisible(set.coresOption(nc))
set.seed(1L)
if (!get.mcOption()) {
  cl <- makeCluster(nc)
  set.ClusterOption(cl)
} else{
  mc.reset.stream()
}
boot_out_par <- aple.mc(as.vector(scale(wheat$yield_detrend, scale=FALSE)),
    nb2listw(nbr12, style="W"), nsim=500)
if (!get.mcOption()) {
  set.ClusterOption(NULL)
  stopCluster(cl)
}
boot_out_par
invisible(set.coresOption(oldCores))
RNGkind(oldRNG[1], oldRNG[2])
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}

