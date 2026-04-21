\name{moran.plot.seismogram}
\alias{moran.plot.seismogram}

\title{Moran seismogram}

\description{
A variant of the Moran scatterplot (see \code{\link{moran.plot}}) supplemented by lines connecting location-wise critical value configurations of attribute values and spatial lags. The plot allows for visual inspection of potential spatial weights misspecifiation.
}
\usage{
moran.plot.seismogram(x, listw, locmoran, alpha = 0.05, adjusted_p = NULL,
  xlab = NULL, ylab = NULL, return_df = TRUE, spChk = NULL,
  zero.policy = attr(listw, "zero.policy"))
}

\arguments{
  \item{x}{a numerical vector holding the attribute of interest}
  \item{listw}{a \code{listw} spatial weights object}
  \item{locmoran}{a fitted object of type localmoran}
  \item{alpha}{default 0.05; the desired significance level regarding local Moran's \emph{I}}
  \item{adjusted_p}{default NULL; an optional vector of \emph{p} values adjusted to account for multiple testing as is returned by \code{\link{p.adjustSP}}; if NULL, standard normal approximation is used to determine critical values}
  \item{xlab}{default NULL; an optional label for the x-axis}
  \item{ylab}{default NULL; an optional label for the y-axis}
  \item{return_df}{default TRUE; invisibly return a data.frame object, if FALSE do not return anything}
  \item{spChk}{default NULL to use \code{get.spChkOption()}; should the locmoran names be checked against the listw spatial objects for identity integrity, TRUE, or FALSE}
  \item{zero.policy}{default option stored in the listw object; if FALSE stop with error for any empty neighbour sets, if TRUE permit the weights list to be formed with zero-length weights vectors}
}

\details{
The Moran seismogram is a version of the Moran scatterplot that is complemented by lines connecting the location-specific critical values of local Moran's \emph{I} in the plot. The y-coordinates associated with the critical value are determined either under the assumption of approximate standard normality of the z-scores of the local Moran's \emph{I} values, or based on \emph{p} values provided through \code{adjusted_p}. In the latter case, the critical value is approximated based on the observation with the highest adjusted \emph{p} value that still satisfies the selected significance level. Lines in quadrants with positive spatial autocorrelation are shown in red and lines in quadrants with negative spatial autocorrelation are shown in blue. The representation is similar to a seismogram for detecting earthquakes and thus reveals potentially suspicious local spatial weights configurations by visualising spikes. The latter are displayed in an integrated manner with their positions in the attribute value range and in connection with the types of the associated spatial patterns (by the quadrants of the scatterplot).
}

\value{
When return_df is TRUE, a data frame object with the following members is returned:
  \item{labels}{either the labels provided or the region ids, if not specified}
  \item{x}{the attribute values}
  \item{wx}{the spatial lags of the centred attribute values}
  \item{b}{the y-coordinates (i.e. hypothetical spatial lags) of the critical values given x}
}

\references{Westerholt, R. (2024): Extending the Moran scatterplot by indications of critical values and \emph{p}-values: introducing the Moran seismogram and the drop plot. Proceedings of the 32nd Annual GIS Research UK Conference (GISRUK), Leeds, UK. \doi{10.5281/zenodo.10897792}}

\author{Rene Westerholt \email{rene.westerholt@tu-dortmund.de}}

\seealso{\code{\link{moran.plot}}}

\examples{
# Boston example (CMEDV; owner-occupied housing in USD)
data(boston)
boston.tr <- sf::st_read(system.file("shapes/boston_tracts.gpkg", package="spData")[1])
boston.nb <- poly2nb(boston.tr)
boston.listw <- nb2listw(boston.nb)
moran.plot.seismogram(boston.c$CMEDV, boston.listw,
 localmoran(boston.c$CMEDV, boston.listw), 0.01, zero.policy = TRUE)
}
