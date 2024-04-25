\name{moran.plot.drop}
\alias{moran.plot.drop}

\title{Moran drop plot}

\description{
A version of the Moran scatterplot (see \code{\link{moran.plot}}) supplemented by lines indicating \emph{p} values for visual inspection of statistical significance.
}
\usage{
moran.plot.drop(x, listw, locmoran, cv = 2.58,
  significant = TRUE, plain = FALSE, zero.policy = FALSE,
  xlab = NULL, ylab = NULL, plot = TRUE, return_df = TRUE,
  spChk = NULL, labels = NULL)
}

\arguments{
  \item{x}{a numerical vector holding the attribute of interest}
  \item{listw}{a \code{listw} spatial weights object}
  \item{locmoran}{a fitted object of type localmoran}
  \item{cv}{default 2.58; the desired critical value assuming approximate normality of standardised local Moran's \emph{I} values}
  \item{significant}{default TRUE; a parameter indicating whether to display plot distances of significant (default) or non-significant observations}
  \item{plain}{default FALSE; a plain Moran scatterplot is displayed if TRUE}
  \item{zero.policy}{default NULL; use global option value; if FALSE stop with error for any empty neighbour sets, if TRUE permit the weights list to be formed with zero-length weights vectors}
  \item{xlab}{label for x axis}
  \item{ylab}{label for y axis}
  \item{plot}{default TRUE; if FALSE, plotting is suppressed}
  \item{return_df}{default TRUE; invisibly return a data.frame object, if FALSE do not return anything}
  \item{spChk}{should the data vector names be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use \code{get.spChkOption()}}
  \item{labels}{default NULL; no labels are plotted by default; character labels for points are assigned to significant observations if provided, region IDs are used if set to TRUE}
}

\details{
The Moran drop plot is a version of the Moran scatterplot supplemented by indications of \emph{p} values. The standard Moran scatterplot provides indirect information about the effect size (distance from the trend line), but not about the level of confidence to determine whether the effects shown might be more than just random outcomes. The Moran drop plot marks significant points in red (positive) and blue (negative spatial autocorrelation) and adds so-called drop lines connecting the significant observations to the scatterplot positions of their associated critical values. The coordinates of the latter are determined under the assumption of an approximate standard normality of the z-scores of the local Moran's \emph{I} values. The longer the lines, the higher the associated \emph{p} values. The visualisation thus enables visual inspection of statistical significance while maintaining the relationship to both attribute value ranges and scatterplot quadrants. It is also possible to invert the visualised relationship and display the distances of non-significant observations to their corresponding critical values (if significant is set to FALSE).
}

\value{
When return_df is TRUE, a data frame object with the following members is returned:
  \item{z}{the standardised attribute values}
  \item{wz}{the standardised spatially lagged attribute values}
  \item{b}{the y-coordinates of the critical values}
  \item{line_lengths}{the absolute distances between b and z}
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
moran.plot.drop(boston.c$CMEDV, boston.listw, localmoran(boston.c$CMEDV, boston.listw), 2.58,
  zero.policy = TRUE, significant = TRUE, plain = FALSE, labels = FALSE, plot = TRUE)
}
