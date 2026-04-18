\name{localGS}
\alias{localGS}
\title{A local hotspot statistic for analysing multiscale datasets}
\description{
The function implements the \eqn{GS_i} test statistic for local hotspots on specific pairwise evaluated distance bands, as proposed by Westerholt et al. (2015).  Like the hotspot estimator \eqn{G_i^*}, the \eqn{GS_i} statistic is given as z-scores that can be evaluated accordingly. The idea of the method is to identify hotspots in datasets that comprise several, difficult-to-separate processes operating at different scales. This is often the case in complex user-generated datasets such as those from Twitter feeds. For example, a football match could be reflected in tweets from pubs, homes, and the stadium vicinity. These exemplified phenomena represent different processes that may be detected at different geometric scales. The \eqn{GS_i} method enables this identification by specifying a geometric scale band and strictly calculating all statistical quantities such as mean and variance solely from respective relevant observations that interact on the range of the adjusted scale band. In addition, in each neighbourhood not only the relationships to the respective central unit, but all scale-relevant relationships are considered. In this way, hotspots can be detected on specific scale ranges independently of other scales. The statistic is given as:
\deqn{GS_i = \frac{\displaystyle\sum_{j; k < j}{w_{ij}w_{ik}\phi_{jk}a_{jk}} - \frac{W_i}{\Phi} \displaystyle\sum_{j; k < j}{\phi_{jk}a_{jk}}}{\sqrt{\frac{W_i}{\Phi}\displaystyle\sum_{j; k < j}{\phi_{jk}a_{jk}^2} + \frac{W_i\left(W_i-1\right)}{\Phi\left(\Phi-1\right)}\left(\Gamma^2 -\!\! \displaystyle\sum_{j; k < j}{\left(\phi_{jk}a_{jk}\right)^2}\right) - \left(\frac{W_i}{\Phi}\displaystyle\sum_{j; k < j}{\phi_{jk}a_{jk}}\right)^2}}}
with
\deqn{a_{jk} = x_j + x_k,\;\;\; W_i = \displaystyle\sum_{j; k < j}{w_{ij}w_{ik}\phi_{jk}},\;\;\; \Phi = \displaystyle\sum_{j; k < j}{\phi_{jk}},\;\;\; \textrm{and} \;\;\; \Gamma = \displaystyle\sum_{j; k < j}{\phi_{jk}a_{jk}}.}

}

\usage{
localGS(x, listw, dmin, dmax, attr, longlat = NULL)
}
\arguments{
  \item{x}{a \code{sf} or \code{sp} object}
  \item{listw}{a \code{listw} object}
  \item{dmin}{a lower distance bound (greater than or equal)}
  \item{dmax}{an upper distance bound (less than or equal)}
  \item{attr}{the name of the attribute of interest}
  \item{longlat}{default NULL; TRUE if point coordinates are longitude-latitude decimal degrees, in which case distances are measured in kilometres; if x is a SpatialPoints object, the value is taken from the object itself, and overrides this argument if not NULL; distances are measured in map units if FALSE or NULL}
}
\details{
Only pairs of observations with a shared distance (in map units) on the interval [\code{dmin}, \code{dmax}] that are within a maximum radius of \code{dmax} around a corresponding output observation are considered. Thereby, also the mean values and variance terms estimated within the measure are adjusted to the scale range under consideration. For application examples of the method see Westerholt et al. (2015) (applied to tweets) and Sonea & Westerholt (2021) (applied in an access to banking scenario).
}
  \value{
    A vector of \eqn{GS_i} values that are given as z-scores.
}
\references{
Westerholt, R., Resch, B. & Zipf, A. 2015. A local scale-sensitive indicator of spatial autocorrelation for assessing high-and low-value clusters in multiscale datasets. International Journal of Geographical Information Science, 29(5), 868--887, \doi{10.1080/13658816.2014.1002499}.

Sonea, A. and Westerholt, R. (2021): Geographic and temporal access to basic banking services in Wales. Applied Spatial Analysis and Policy, 14 (4), 879--905, \doi{10.1007/s12061-021-09386-3}.
}

\author{René Westerholt \email{rene.westerholt@tu-dortmund.de}}

\seealso{\code{\link{localG}}}
  
\examples{
\donttest{
boston.tr <- sf::st_read(system.file("shapes/boston_tracts.gpkg", package="spData")[1])
boston.tr_utm <- st_transform(boston.tr, 32619) #26786

boston_listw1 <- nb2listwdist(dnearneigh(st_centroid(boston.tr_utm), 1, 2000),
    boston.tr_utm, type = "dpd", alpha = 2, zero.policy = TRUE, dmax = 9500)

boston_listw2 <- nb2listwdist(dnearneigh(st_centroid(boston.tr_utm), 5000, 9500), 
    boston.tr_utm, type = "dpd", alpha = 2, zero.policy = TRUE, dmax = 9500)

boston_RM_gsi_1 <- localGS(boston.tr_utm, boston_listw1, 1, 2000, "RM", FALSE)
boston_RM_gsi_2 <- localGS(boston.tr_utm, boston_listw2, 2000, 9500, "RM", FALSE)
}
}
\keyword{spatial}
