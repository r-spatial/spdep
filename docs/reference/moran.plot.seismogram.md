# Moran seismogram

A variant of the Moran scatterplot (see
[`moran.plot`](https://r-spatial.github.io/spdep/reference/moran.plot.md))
supplemented by lines connecting location-wise critical value
configurations of attribute values and spatial lags. The plot allows for
visual inspection of potential spatial weights misspecifiation.

## Usage

``` r
moran.plot.seismogram(x, listw, locmoran, alpha = 0.05, adjusted_p = NULL,
  xlab = NULL, ylab = NULL, return_df = TRUE, spChk = NULL,
  zero.policy = attr(listw, "zero.policy"), na.action = na.fail,
  conditional = TRUE, alternative = "two.sided", mlvar = TRUE,
  adjust.x = FALSE, quadrant.type = "mean")
```

## Arguments

- x:

  a numerical vector holding the attribute of interest

- listw:

  a `listw` spatial weights object

- locmoran:

  (may be omitted, computed internally) a fitted object of class
  localmoran; if NULL, computed internally

- alpha:

  default 0.05; the desired two-sided significance level regarding local
  Moran's *I*

- adjusted_p:

  default NULL; an optional vector of *p* values adjusted to account for
  multiple testing as is returned by
  [`p.adjustSP`](https://r-spatial.github.io/spdep/reference/p.adjustSP.md);
  if NULL, standard normal approximation is used to determine critical
  values

- xlab:

  default NULL; an optional label for the x-axis

- ylab:

  default NULL; an optional label for the y-axis

- return_df:

  default TRUE; invisibly return a data.frame object, if FALSE do not
  return anything

- spChk:

  default NULL to use
  [`get.spChkOption()`](https://r-spatial.github.io/spdep/reference/set.spChkOption.md);
  should the locmoran names be checked against the listw spatial objects
  for identity integrity, TRUE, or FALSE

- zero.policy:

  default option stored in the listw object; if FALSE stop with error
  for any empty neighbour sets, if TRUE permit the weights list to be
  formed with zero-length weights vectors

- na.action:

  a function (default `na.fail`), can also be `na.omit` or
  `na.exclude` - in these cases the weights list will be subsetted to
  remove NAs in the data. It may be necessary to set zero.policy to TRUE
  because this subsetting may create no-neighbour observations. Note
  that only weights lists created without using the glist argument to
  `nb2listw` may be subsetted. If `na.pass` is used, zero is substituted
  for NA values in calculating the spatial lag. (Note that na.exclude
  will only work properly starting from R 1.9.0, na.omit and na.exclude
  assign the wrong classes in 1.8.\*)

- conditional:

  default TRUE: expectation and variance are calculated using the
  conditional randomization null (Sokal 1998 Eqs. A7 & A8). Elaboration
  of these changes available in Sauer et al. (2021). If FALSE:
  expectation and variance are calculated using the total randomization
  null (Sokal 1998 Eqs. A3 & A4).

- alternative:

  a character string specifying the alternative hypothesis, must be one
  of greater, less or two.sided (default).

- mlvar:

  default TRUE: values of local Moran's I are reported using the
  variance of the variable of interest (sum of squared deviances over
  n), but can be reported as the sample variance, dividing by (n-1)
  instead; both are used in other implementations.

- adjust.x:

  default FALSE, if TRUE, x values of observations with no neighbours
  are omitted in the mean of x

- quadrant.type:

  Default `"mean"`, for `"localmoran"` objects only, can be
  `c("mean", "median", "pysal")` to partition the Moran scatterplot;
  `"mean"` partitions on the means of the variable and its spatial lag,
  `"median"` on medians of the variable and its spatial lag, `"pysal"`
  at zero for the centred variable and its spatial lag

## Details

The Moran seismogram is a version of the Moran scatterplot that is
complemented by lines connecting the location-specific critical values
of local Moran's *I* in the plot. The y-coordinates associated with the
critical value are determined either under the assumption of approximate
standard normality of the z-scores of the local Moran's *I* values, or
based on *p* values provided through `adjusted_p`. In the latter case,
the critical value is approximated based on the observation with the
highest adjusted *p* value that still satisfies the selected
significance level. Lines in quadrants with positive spatial
autocorrelation are shown in red and lines in quadrants with negative
spatial autocorrelation are shown in blue. The representation is similar
to a seismogram for detecting earthquakes and thus reveals potentially
suspicious local spatial weights configurations by visualising spikes.
The latter are displayed in an integrated manner with their positions in
the attribute value range and in connection with the types of the
associated spatial patterns (by the quadrants of the scatterplot).

## Value

When return_df is TRUE, a data frame object with the following members
is returned:

- labels:

  either the labels provided or the region ids, if not specified

- x:

  the attribute values

- wx:

  the spatial lags of the centred attribute values

- b:

  the y-coordinates (i.e. hypothetical spatial lags) of the critical
  values given x

## References

Westerholt, R. (2024): Extending the Moran scatterplot by indications of
critical values and *p*-values: introducing the Moran seismogram and the
drop plot. Proceedings of the 32nd Annual GIS Research UK Conference
(GISRUK), Leeds, UK.
[doi:10.5281/zenodo.10897792](https://doi.org/10.5281/zenodo.10897792)

## Author

Rene Westerholt <rene.westerholt@tu-dortmund.de>

## See also

[`moran.plot`](https://r-spatial.github.io/spdep/reference/moran.plot.md),
[`localmoran`](https://r-spatial.github.io/spdep/reference/localmoran.md),
[`hotspot.localmoran`](https://r-spatial.github.io/spdep/reference/hotspotmap.md)

## Examples

``` r
# Boston example (CMEDV; owner-occupied housing in USD)
data(boston)
boston.tr <- sf::st_read(system.file("shapes/boston_tracts.gpkg", package="spData")[1])
#> Reading layer `boston_tracts' from data source 
#>   `/home/rsb/lib/r_libs/spData/shapes/boston_tracts.gpkg' using driver `GPKG'
#> Simple feature collection with 506 features and 36 fields
#> Geometry type: POLYGON
#> Dimension:     XY
#> Bounding box:  xmin: -71.52311 ymin: 42.00305 xmax: -70.63823 ymax: 42.67307
#> Geodetic CRS:  NAD27
boston.nb <- poly2nb(boston.tr)
boston.listw <- nb2listw(boston.nb)
moran.plot.seismogram(boston.c$CMEDV, boston.listw, alpha = 0.01,
 zero.policy = TRUE)
```
