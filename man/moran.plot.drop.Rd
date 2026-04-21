\name{moran.plot.drop}
\alias{moran.plot.drop}

\title{Moran drop plot}

\description{
A version of the Moran scatterplot (see \code{\link{moran.plot}}) supplemented by lines indicating \emph{p} values for visual inspection of statistical significance.
}
\usage{
moran.plot.drop(x, listw, locmoran, alpha = 0.05, adjusted_p = NULL,
  significant = TRUE, xlab = NULL, ylab = NULL, return_df = TRUE,
  spChk = NULL, labels = NULL, zero.policy=attr(listw, "zero.policy"))
}

\arguments{
  \item{x}{a numerical vector holding the attribute of interest}
  \item{listw}{a \code{listw} spatial weights object}
  \item{locmoran}{a fitted object of type \code{localmoran}}
  \item{alpha}{default 0.05; the desired significance level regarding local Moran's \emph{I}}
  \item{adjusted_p}{default NULL; an optional vector of \emph{p} values adjusted to account for multiple testing as is returned by \code{\link{p.adjustSP}}; if NULL, standard normal approximation is used to determine critical values}
  \item{significant}{default TRUE; a parameter indicating whether to display critical value distances of significant (default) or non-significant observations}
  \item{xlab}{default NULL; an optional label for the x-axis}
  \item{ylab}{default NULL; an optional label for the y-axis}
  \item{return_df}{default TRUE; invisibly return a data.frame object, if FALSE do not return anything}
  \item{spChk}{default NULL to use \code{get.spChkOption()}; should the locmoran names be checked against the listw spatial objects for identity integrity, TRUE, or FALSE}
  \item{labels}{default NULL; no labels are plotted by default; region IDs from the \code{listw} object are used as labels of significant observations if set to TRUE; custom labels are used if a character vector is provided}
  \item{zero.policy}{default option stored in the listw object; if FALSE stop with error for any empty neighbour sets, if TRUE permit the weights list to be formed with zero-length weights vectors}
}

\details{
The Moran drop plot is a version of the Moran scatterplot supplemented by visual indications of \emph{p} values. The standard Moran scatterplot provides information about the effect size but not about the level of confidence to determine whether the effects shown might be more than just random outcomes. The Moran drop plot marks significant points in red (positive) and blue (negative spatial autocorrelation) and adds so-called drop lines connecting the significant observations to the scatterplot positions of their associated critical values. The coordinates of the latter are determined either under the assumption of approximate standard normality of the z-scores of the local Moran's \emph{I} values, or based on \emph{p} values provided through \code{adjusted_p}. In the latter case, the critical value is approximated based on the observation with the highest adjusted \emph{p} value that still satisfies the selected significance level. The longer the lines, the lower the associated \emph{p} values. The visualisation thus enables visual inspection of statistical significance while maintaining the relationship to both attribute value ranges and scatterplot quadrants. It is also possible to invert the visualised relationship and display the distances of non-significant observations to their corresponding critical values (if significant is set to FALSE).
}

\value{
When return_df is TRUE, a data frame object with the following members is returned:
  \item{labels}{either the labels provided or the region ids, if not specified}
  \item{x}{the attribute values}
  \item{wx}{the spatial lags of the centred attribute values}
  \item{b}{the y-coordinates (i.e. hypothetical spatial lags) of the critical values given x}
  \item{line_lengths}{the absolute distances between b and x}
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
moran.plot.drop(boston.c$CMEDV, boston.listw, localmoran(boston.c$CMEDV, boston.listw), 0.01,
  significant = TRUE, labels = NULL)
}
