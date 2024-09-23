\name{spatialdelta}
\alias{spatialdelta}
\alias{summary.spatialdelta}
\alias{print.summary.spatialdelta}
\alias{linearised_diffusive_weights}
\alias{metropolis_hastings_weights}
\alias{iterative_proportional_fitting_weights}
\alias{graph_distance_weights}
\alias{plot_spatialcoords}
\alias{plot_spatialcoords.default}
\alias{plot_spatialcoords.spatialdelta}
\alias{plot_moran}
\alias{plot_moran.default}
\alias{plot_moran.spatialdelta}
\alias{plot_spatialscree}
\alias{plot_spatialscree.default}
\alias{plot_spatialscree.spatialdelta}
\alias{factorial_coordinates}
\alias{factorial_coordinates.default}
\alias{factorial_coordinates.default}
\alias{factorial_coordinates.spatialdelta}
\alias{plot_factorialcoords}
\alias{plot_factorialcoords.default}
\alias{plot_factorialcoords.spatialdelta}
\alias{plot_factorialscree}
\alias{plot_factorialscree.default}
\alias{plot_factorialscree.spatialdelta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Weighted Multivariate Spatial Autocorrelation Measures}
\description{The kernel-based weighted multivariate spatial autocorrelation measure delta proposed in Bavaud (2024) is implemented, together with support functions to create spatial weights matrices from symmetric binary adjacency matrices; summary and display methods are provided.}
\usage{
spatialdelta(D, w, f, alternative = "greater", CF = FALSE)
\method{summary}{spatialdelta}(object, ...)
\method{print}{summary.spatialdelta}(x, digits=getOption("digits"), ...)
linearised_diffusive_weights(A, f, t_choice = 2)
metropolis_hastings_weights(A, f)
iterative_proportional_fitting_weights(A, f, g=0.001, iter=1000, tol=1e-10,
 tol.margins=1e-10, print=FALSE)
graph_distance_weights(A, f)
plot_spatialcoords(x, ...)
\method{plot_spatialcoords}{default}(x, ...)
\method{plot_spatialcoords}{spatialdelta}(x, cols = c(1L, 2L),
 mult = c(1, 1), power = 1L, fmult = NULL, names = NULL, bg = 1,
 pos = 3, cex = 0.6, ...)
plot_moran(x, y, ...)
\method{plot_moran}{default}(x, y, ...)
\method{plot_moran}{spatialdelta}(x, y, fmult = NULL, names = NULL, bg = 1,
 pos = 3, cex = 0.6, ...)
plot_spatialscree(x, ...)
\method{plot_spatialscree}{default}(x, ...)
\method{plot_spatialscree}{spatialdelta}(x, ...)
factorial_coordinates(x)
\method{factorial_coordinates}{default}(x)
\method{factorial_coordinates}{spatialdelta}(x)
plot_factorialcoords(x, ...)
\method{plot_factorialcoords}{default}(x, ...)
\method{plot_factorialcoords}{spatialdelta}(x, cols = c(1L, 2L),
 mult = c(1, 1), fmult = NULL, names = NULL, bg = 1, pos = 3, cex = 0.6, ...)
plot_factorialscree(x, ...)
\method{plot_factorialscree}{default}(x, ...)
\method{plot_factorialscree}{spatialdelta}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{D}{\code{D} is a square matrix of dissimilarities between (possibly multivariate) observations}
  \item{w}{\code{w} is a square matrix of spatial weights as returned by construction functions \code{linearised_diffusive_weights}, \code{metropolis_hastings_weights}, \code{iterative_proportional_fitting_weights}, \code{graph_distance_weights} or similar}
  \item{f}{\code{f} are weights reflecting the contribution of each observation to the (possibly multivariate) data set, they may be uniform}
  \item{alternative}{a character string specifying the alternative hypothesis, must be one of "greater" (default), "less" or "two.sided"}
  \item{CF}{default \code{FALSE}; if \code{TRUE}, use a Cornish-Fisher correction to the standard deviate}
  \item{object, x}{\code{object}, \code{x} are objects returned by \code{spatialdelta} of class \code{"spatialdelta"}}
  \item{digits}{default \code{getOption("digits")}, or a non-null value specifying the minimum number of significant digits to be printed in values}
  \item{A}{\code{A} is a symmetric binary adjacency matrix}
  \item{t_choice}{default 2, the inverse of the largest eigenvalue of the adjusted Laplacian adjacency matrix (t2, Bavaud 2024, page 583), otherwise 1 (t1, page 579)}
  \item{g}{default 0.001, (Bavaud 2024, page 579, 589), a small quantity to lift binary adjacency matrix above zero}
  \item{iter, tol, tol.margins, print}{arguments passed to \code{\link[mipfp]{Ipfp}}}
  \item{y}{a numeric vector from which to create a Moran plot}
  \item{cols}{integer vector of length 2, default \code{c(1L, 2L)}, giving the two columns of the regional or factorial coordinates to plot}
  \item{mult}{numeric vector of length 2, default \code{c(1, 1)}, to scale or reverse axes as signs are not given to be the same for different eigenproblem implementations}
  \item{power}{integer vector of length 1, default 1; if higher, use a powered kernel}
  \item{fmult}{default \code{NULL} for automatic scaling of the circle symbols expressing the observation weights \code{f} in the plot as \code{(0.02*diff(range(X)))/diff(range(f))} where \code{X} is the first vector of coordinates; may be set manually to a numeric value}
  \item{names}{default \code{NULL} yielding integer index values, if a character vector is given, the coordinates will be annotated with those strings rather than indices}
  \item{bg, pos, cex}{arguments passed to graphics functions or methods}
  \item{...}{other arguments passed to other functions or methods}
}
\details{In general, the user should choose one of the functions for constructing spatial weights from the symmetric binary, or similar to symmetric - for example row-standardised - adjacency matrix.}

\value{\code{spatialdelta} returns an object with classes \code{"spatialdelta"} and \code{"htest"}. Its \code{print} method treats it as an \code{"htest"} object. Other methods use additional attributes:
  \item{Vd0}{Variance of delta calculated in another way}
  \item{Kx}{Multivariate kernel matrix}
  \item{Kw}{Spatial weights kernel matrix}
  \item{f}{Input observation weights}
  \item{w}{Input spatial weights}
  \item{VI}{Spatial component of delta variance}
  \item{Vx}{Multivariate component of delta variance}
  \item{alphamu}{Spatial component of delta skewness}
  \item{alphalambda}{Multivariate component of delta skewness}
  \item{gammamu}{Spatial component of delta excess kurtosis}
  \item{gammalambda}{Multivariate component of delta excess kurtosis}
\code{linearised_diffusive_weights}, \code{metropolis_hastings_weights}, \code{iterative_proportional_fitting_weights} and \code{graph_distance_weights} return square spatial weights matrices (the latter two return matrices created by \code{metropolis_hastings_weights} if \pkg{igraph} or \pkg{mipfp} are not available); \code{factorial_coordinates} returns a square factorial matrix. 
}
\references{Bavaud, F. (2024), Measuring and Testing Multivariate Spatial Autocorrelation in a Weighted Setting: A Kernel Approach. Geographical Analysis, 56: 573-599. \doi{10.1111/gean.12390}}
\author{Roger Bivand, \email{Roger.Bivand@nhh.no}}
\note{The input adjacency matrix should be symmetric, and not have sub-graphs or islands - it should be possible to step across graph edges from any observation to any other.}

\seealso{\code{\link[mipfp]{Ipfp}}, \code{\link[igraph]{distances}}}
\examples{
toy_f <- c(0.4, 0.1, 0.3, 0.2)
toy_xa <- c(4, 9, 1, 1)
toy_nba <- structure(list(c(2L, 3L, 4L), c(1L), c(1L, 4L), c(1L, 3L)), class="nb")
toy_gla <- list(c(0.25, 0.5, 0.25), c(1), c(2/3, 1/3), c(0.5, 0.5))
toy_wa <- spdep::nb2mat(toy_nba, glist=toy_gla, style="W")
spatialdelta(as.matrix(dist(toy_xa))^2, toy_wa, toy_f, alternative="two.sided")
toy_xb <- toy_xa
toy_glb <- list(c(0.5, 0.125, 0.25, 0.125), c(0.5, 0.5), c(1/3, 0.5, 1/6), c(0.25, 0.25, 0.5))
toy_nbb <- spdep::include.self(toy_nba)
toy_wb <- spdep::nb2mat(toy_nbb, glist=toy_glb, style="W")
spatialdelta(as.matrix(dist(toy_xb))^2, toy_wb, toy_f, alternative="two.sided")
toy_xc <- c(3, 3, 1, 1)
toy_nbc <- structure(list(c(3L, 4L), c(3L, 4L), c(1L, 2L), c(1L, 2L)), class="nb")
toy_glc <- list(c(0.6, 0.4), c(0.6, 0.4), c(0.8, 0.2), c(0.8, 0.2))
toy_wc <- spdep::nb2mat(toy_nbc, glist=toy_glc, style="W")
spatialdelta(as.matrix(dist(toy_xc))^2, toy_wc, toy_f, alternative="two.sided")
\donttest{
library(Guerry)
data(gfrance85)
gfrance85 <- st_as_sf(gfrance85)
moral <- scale(st_drop_geometry(gfrance85)[, 7:12])
Dmoral <- as.matrix(dist(moral))^2
fG <- gfrance85$Pop1831/sum(gfrance85$Pop1831)
gnb <- poly2nb(gfrance85, row.names=gfrance85$Department)
frG <- nb2mat(gnb, style="B")
ldwG <- linearised_diffusive_weights(frG, fG)
mhwG <- metropolis_hastings_weights(frG, fG)
ifwG <- iterative_proportional_fitting_weights(frG, fG)
gdwG <- graph_distance_weights(frG, fG)
moral_ldw <- spatialdelta(Dmoral, ldwG, fG, CF=FALSE)
moral_mhw <- spatialdelta(Dmoral, mhwG, fG, CF=FALSE)
moral_ifw <- spatialdelta(Dmoral, ifwG, fG, CF=FALSE)
moral_gdw <- spatialdelta(Dmoral, gdwG, fG, CF=FALSE)
opar <- par(mfrow=c(1, 4))
plot_spatialcoords(moral_ldw, cols=c(2,1), mult=c(-1, -1),
 names=gfrance85$Department, main="Linearised Diffusive", asp=1)
plot_spatialcoords(moral_mhw, cols=c(2,1), mult=c(1, 1),
 names=gfrance85$Department, main="Metropolis-Hastings", asp=1)
plot_spatialcoords(moral_ifw, cols=c(2,1), mult=c(1, 1),
 names=gfrance85$Department, main="Iterative fitting", asp=1)
plot_spatialcoords(moral_gdw, cols=c(2,1), mult=c(1, 1),
 names=gfrance85$Department, main="Graph distance", asp=1)
par(opar)
opar <- par(mfrow=c(1, 4))
plot_spatialscree(moral_ldw, main="linearised diffusive")
plot_spatialscree(moral_mhw, main="Metropolis-Hastings")
plot_spatialscree(moral_ifw, main="iterative fitting")
plot_spatialscree(moral_gdw, main="graph distance")
par(opar)
Crime_pers <- scale(gfrance85$Crime_pers)
funif <- rep(1/length(Crime_pers), length(Crime_pers))
lw <- nb2listw(gnb, style="W")
spatialdelta(as.matrix(dist(Crime_pers))^2, listw2mat(lw), funif)$estimate["delta"]
moran.test(Crime_pers, lw, randomisation=FALSE)$estimate["Moran I statistic"]
}
}
\keyword{spatial}
\keyword{multivariate}

