\name{moran.plot.seismogram}
\alias{moran.plot.seismogram}

\title{Moran seismogram}

\description{
A variant of the Moran scatterplot (see \code{\link{moran.plot}}) supplemented by lines connecting location-wise critical values. The plot allows for visual inspection of potential spatial weights misspecifiation.
}
\usage{
moran.plot.seismogram(x, listw, locmoran, cv = 2.58,
  plain = FALSE, zero.policy = FALSE, xlab = NULL,
  ylab = NULL, plot = TRUE, return_df = TRUE, spChk = NULL)
}

\arguments{
  \item{x}{a numerical vector holding the attribute of interest}
  \item{listw}{a \code{listw} spatial weights object}
  \item{locmoran}{a fitted object of type localmoran}
  \item{cv}{default 2.58; the desired critical value assuming approximate normality of standardised local Moran's \emph{I} values}
  \item{plain}{default FALSE; a plain Moran scatterplot is displayed if TRUE}
  \item{zero.policy}{default NULL; use global option value; if FALSE stop with error for any empty neighbour sets, if TRUE permit the weights list to be formed with zero-length weights vectors}
  \item{xlab}{label for x axis}
  \item{ylab}{label for y axis}
  \item{plot}{default TRUE; if FALSE, plotting is suppressed}
  \item{return_df}{default TRUE; invisibly return a data.frame object, if FALSE do not return anything}
  \item{spChk}{should the data vector names be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use \code{get.spChkOption()}}
}

\details{
The Moran seismogram is a version of the Moran scatterplot that is complemented by lines connecting the location-specific critical values of local Moran's \emph{I}. Lines in quadrants with positive spatial autocorrelation are shown in red and lines in quadrants with negative spatial autocorrelation are shown in blue. The representation is similar to a seismogram for detecting earthquakes and thus reveals potentially suspicious local spatial weights configurations by visualising spikes. The latter are displayed in an integrated manner with their positions in the attribute value range and in connection with the types of the associated spatial patterns (by the quadrants of the scatterplot).
}

\value{
When return_df is TRUE, a data frame object with the following members is returned:
  \item{z}{the standardised attribute values}
  \item{wz}{the standardised spatially lagged attribute values}
  \item{b}{the y-coordinates of the critical values}
}

\references{Westerholt, R. (2024): Extending the Moran scatterplot by indications of critical values and \emph{p}-values: introducing the Moran seismogram and the drop plot. Proceedings of the 32nd Annual GIS Research UK Conference (GISRUK), Leeds, UK. \url{https://doi.org/10.5281/zenodo.10897792}}

\author{Rene Westerholt \email{rene.westerholt@tu-dortmund.de}}

\seealso{\code{\link{moran.plot}}}

\examples{
# Boston example (CMEDV; owner-occupied housing in USD)
data(boston)
boston.tr <- sf::st_read(system.file("shapes/boston_tracts.shp", package="spData")[1])
boston.nb <- poly2nb(boston.tr)
boston.listw <- nb2listw(boston.nb)
moran.plot.seismogram(boston.c$CMEDV, boston.listw, localmoran(boston.c$CMEDV, boston.listw), 2.58,
zero.policy = TRUE, plain = FALSE, plot = TRUE)
}
