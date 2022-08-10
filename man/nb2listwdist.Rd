% Copyright 2020 by Roger S. Bivand and Rene Westerholt
\name{nb2listwdist}
\alias{nb2listwdist}
\title{Distance-based spatial weights for neighbours lists}
\description{
The \code{nb2listwdist} function supplements a neighbours list with spatial weights for the chosen types of distance modelling and coding scheme. While the offered coding schemes parallel those of the \code{nb2listw} function, three distance-based types of weights are available: inverse distance weighting (IDW), double-power distance weights, and exponential distance decay. The \code{can.be.simmed} helper function checks whether a spatial weights object is similar to symmetric and can be so transformed to yield real eigenvalues or for Cholesky decomposition.
}
\usage{
nb2listwdist(neighbours, x, type="idw", style="raw", 
  alpha = 1, dmax = NULL, longlat = NULL, zero.policy=NULL)
}
\arguments{
  \item{neighbours}{an object of class \code{nb}}
  \item{x}{an \code{sp} \code{sf}, or \code{sfc} object}
  \item{type}{default \dQuote{idw}; the intended type of distance modelling, can take values \dQuote{idw}, \dQuote{exp}, and \dQuote{dpd}}
  \item{style}{default \dQuote{raw}; \code{style} can take values \dQuote{raw}, \dQuote{W}, \dQuote{B}, \dQuote{C}, \dQuote{U}, \dQuote{minmax}, and \dQuote{S}}
  \item{alpha}{default 0; a parameter for controlling the distance modelling, see \dQuote{Details}}
  \item{dmax}{default NULL, maximum distance threshold that is required for type \dQuote{dpd} but optional for all other types}
  \item{longlat}{default NULL; TRUE if point coordinates are longitude-latitude decimal degrees, in which case distances are measured in metres; if x is a SpatialPoints object, the value is taken from the object itself, and overrides this argument if not NULL; distances are measured in map units if FALSE or NULL}
  \item{zero.policy}{default NULL; use global option value; if FALSE stop with error for any empty neighbour sets, if TRUE permit the weights list to be formed with zero-length weights vectors}
}
\details{
Starting from a binary neighbours list, in which regions are either listed as neighbours or are absent (thus not in the set of neighbours for some definition), the function adds a distance-based weights list. Three types of distance weight calculations based on pairwise distances \eqn{d_{ij}} are possible, all of which are controlled by parameter \dQuote{alpha} (\eqn{\alpha} below):
\deqn{\textrm{idw: } w_{ij} = d_{ij}^{-\alpha},}
\deqn{\textrm{exp: } w_{ij} = \exp(-\alpha \cdot d_{ij}),}
\deqn{\textrm{dpd: } w_{ij} = \left[1 - \left(d_{ij}/d_{\textrm{max}}\right)^{\alpha}\right]^{\alpha},}
the latter of which leads to \eqn{w_{ij} = 0} for all \eqn{d_{ij} > d_{\textrm{max}}}. Note that IDW weights show extreme behaviour close to 0 and can take on the value infinity. In such cases, the infinite values are replaced by the largest finite weight present in the weights list.

The default coding scheme is \dQuote{raw}, which outputs the raw distance-based weights without applying any kind of normalisation. In addition, the same coding scheme styles that are also available in the \code{nb2listw} function can be chosen. B is the basic binary coding, W is row standardised (sums over all links to n), C is globally standardised (sums over all links to n), U is equal to C divided by the number of neighbours (sums over all links to unity), while S is the variance-stabilising coding scheme proposed by Tiefelsdorf et al. 1999, p. 167-168 (sums over all links to n). The \dQuote{minmax} style is based on Kelejian and Prucha (2010), and divides the weights by the minimum of the maximum row sums and maximum column sums of the input weights. It is similar to the C and U styles; it is also available in Stata.

If zero.policy is set to TRUE, weights vectors of zero length are inserted for regions without neighbour in the neighbours list. These will in turn generate lag values of zero, equivalent to the sum of products of the zero row \code{t(rep(0, length=length(neighbours))) \%*\% x}, for arbitraty numerical vector \code{x} of length \code{length(neighbours)}. The spatially lagged value of x for the zero-neighbour region will then be zero, which may (or may not) be a sensible choice.
}                                                                                                  
\value{
A \code{listw} object with the following members:
  \item{style}{one of W, B, C, U, S, minmax as above}
  \item{type}{one of idw, exp, dpd as above}
  \item{neighbours}{the input neighbours list}
  \item{weights}{the weights for the neighbours and chosen style, with attributes set to report the type of relationships (binary or general, if general the form of the glist argument), and style as above}
}
\references{Tiefelsdorf, M., Griffith, D. A., Boots, B. 1999 A variance-stabilizing coding scheme for spatial link matrices, Environment and Planning A, 31, pp. 165--180; Kelejian, H. H., and I. R. Prucha. 2010. Specification and estimation of spatial autoregressive models with autoregressive and heteroskedastic disturbances. Journal of Econometrics, 157: pp. 53--67.}

\author{Rene Westerholt \email{rene.westerholt@tu-dortmund.de}}

\seealso{\code{\link{nb2listw}}, \code{\link{summary.nb}}}

\examples{
# World examples
data(world, package="spData")
# neighbours on distance interval [0, 1000] kilometres
# suppressWarnings(st_crs(world) <- "+proj=longlat") # for older PROJ
pts <- st_centroid(st_transform(world, 3857))
nb_world <- dnearneigh(pts, 0, 1000000)
# Moran's I (life expectancy) with IDW with alpha = 2, no coding scheme
world_weights <- nb2listwdist(nb_world, as(pts, "Spatial"), type = "idw",
  alpha = 2, zero.policy = TRUE)
moran.test(world$lifeExp, world_weights, zero.policy = TRUE, na.action = na.pass)
\dontrun{
# Moran's I (life expectancy) with IDW with alpha = 2, no coding scheme
world_weights <- nb2listwdist(nb_world, pts, type = "idw",
  alpha = 2, zero.policy = TRUE)
moran.test(world$lifeExp, world_weights, zero.policy = TRUE, na.action = na.pass)
# Moran's I (life expectancy), DPD, alpha = 2, dmax = 1000 km, no coding scheme
world_weights <- nb2listwdist(nb_world, pts, type = "dpd",
  dmax = 1000000, alpha = 2, zero.policy = TRUE)
moran.test(world$lifeExp, world_weights, zero.policy = TRUE, na.action = na.pass)
# Boston examples
data(boston, package="spData")
boston_coords <- data.frame(x = boston.utm[,1], y = boston.utm[,2])
boston.geoms <- st_as_sf(boston_coords, coords = c("x", "y"), remove = FALSE)
nb_boston <- dnearneigh(boston.geoms, 0, 3)
# Moran's I (crime) with exp weights with alpha = 2, no coding scheme
boston_weights <- nb2listwdist(nb_boston, boston.geoms, type = "exp", alpha = 2,
  style="raw", zero.policy = TRUE)
moran.test(boston.c$CRIM, boston_weights, zero.policy = TRUE, na.action = na.pass)
# Moran's I (crime) with idw weights with alpha = 2, coding scheme = W
boston_weights <- nb2listwdist(nb_boston, boston.geoms, type = "idw", alpha = 2,
  style="W", zero.policy = TRUE)
moran.test(boston.c$CRIM, boston_weights, zero.policy = TRUE, na.action = na.pass)
}
}

\keyword{spatial}
