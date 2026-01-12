# A local hotspot statistic for analysing multiscale datasets

The function implements the \\GS_i\\ test statistic for local hotspots
on specific pairwise evaluated distance bands, as proposed by Westerholt
et al. (2015). Like the hotspot estimator \\G_i^\*\\, the \\GS_i\\
statistic is given as z-scores that can be evaluated accordingly. The
idea of the method is to identify hotspots in datasets that comprise
several, difficult-to-separate processes operating at different scales.
This is often the case in complex user-generated datasets such as those
from Twitter feeds. For example, a football match could be reflected in
tweets from pubs, homes, and the stadium vicinity. These exemplified
phenomena represent different processes that may be detected at
different geometric scales. The \\GS_i\\ method enables this
identification by specifying a geometric scale band and strictly
calculating all statistical quantities such as mean and variance solely
from respective relevant observations that interact on the range of the
adjusted scale band. In addition, in each neighbourhood not only the
relationships to the respective central unit, but all scale-relevant
relationships are considered. In this way, hotspots can be detected on
specific scale ranges independently of other scales. The statistic is
given as: \$\$GS_i = \frac{\displaystyle\sum\_{j; k \<
j}{w\_{ij}w\_{ik}\phi\_{jk}a\_{jk}} - \frac{W_i}{\Phi}
\displaystyle\sum\_{j; k \<
j}{\phi\_{jk}a\_{jk}}}{\sqrt{\frac{W_i}{\Phi}\displaystyle\sum\_{j; k \<
j}{\phi\_{jk}a\_{jk}^2} +
\frac{W_i\left(W_i-1\right)}{\Phi\left(\Phi-1\right)}\left(\Gamma^2
-\\\\ \displaystyle\sum\_{j; k \<
j}{\left(\phi\_{jk}a\_{jk}\right)^2}\right) -
\left(\frac{W_i}{\Phi}\displaystyle\sum\_{j; k \<
j}{\phi\_{jk}a\_{jk}}\right)^2}}\$\$ with \$\$a\_{jk} = x_j + x_k,\\\\\\
W_i = \displaystyle\sum\_{j; k \< j}{w\_{ij}w\_{ik}\phi\_{jk}},\\\\\\
\Phi = \displaystyle\sum\_{j; k \< j}{\phi\_{jk}},\\\\\\ \textrm{and}
\\\\\\ \Gamma = \displaystyle\sum\_{j; k \< j}{\phi\_{jk}a\_{jk}}.\$\$

## Usage

``` r
localGS(x, listw, dmin, dmax, attr, longlat = NULL)
```

## Arguments

- x:

  a `sf` or `sp` object

- listw:

  a `listw` object

- dmin:

  a lower distance bound (greater than or equal)

- dmax:

  an upper distance bound (less than or equal)

- attr:

  the name of the attribute of interest

- longlat:

  default NULL; TRUE if point coordinates are longitude-latitude decimal
  degrees, in which case distances are measured in kilometres; if x is a
  SpatialPoints object, the value is taken from the object itself, and
  overrides this argument if not NULL; distances are measured in map
  units if FALSE or NULL

## Details

Only pairs of observations with a shared distance (in map units) on the
interval \[`dmin`, `dmax`\] that are within a maximum radius of `dmax`
around a corresponding output observation are considered. Thereby, also
the mean values and variance terms estimated within the measure are
adjusted to the scale range under consideration. For application
examples of the method see Westerholt et al. (2015) (applied to tweets)
and Sonea & Westerholt (2021) (applied in an access to banking
scenario).

## Value

A vector of \\GS_i\\ values that are given as z-scores.

## References

Westerholt, R., Resch, B. & Zipf, A. 2015. A local scale-sensitive
indicator of spatial autocorrelation for assessing high-and low-value
clusters in multiscale datasets. International Journal of Geographical
Information Science, 29(5), 868–887,
[doi:10.1080/13658816.2014.1002499](https://doi.org/10.1080/13658816.2014.1002499)
.

Sonea, A. and Westerholt, R. (2021): Geographic and temporal access to
basic banking services in Wales. Applied Spatial Analysis and Policy, 14
(4), 879–905,
[doi:10.1007/s12061-021-09386-3](https://doi.org/10.1007/s12061-021-09386-3)
.

## Author

René Westerholt <rene.westerholt@tu-dortmund.de>

## See also

[`localG`](https://r-spatial.github.io/spdep/reference/localG.md)

## Examples

``` r
# \donttest{
boston.tr <- sf::st_read(system.file("shapes/boston_tracts.gpkg", package="spData")[1])
#> Reading layer `boston_tracts' from data source 
#>   `/home/rsb/lib/r_libs/spData/shapes/boston_tracts.gpkg' using driver `GPKG'
#> Simple feature collection with 506 features and 36 fields
#> Geometry type: POLYGON
#> Dimension:     XY
#> Bounding box:  xmin: -71.52311 ymin: 42.00305 xmax: -70.63823 ymax: 42.67307
#> Geodetic CRS:  NAD27
boston.tr_utm <- st_transform(boston.tr, 32619) #26786

boston_listw1 <- nb2listwdist(dnearneigh(st_centroid(boston.tr_utm), 1, 2000),
    boston.tr_utm, type = "dpd", alpha = 2, zero.policy = TRUE, dmax = 9500)
#> Warning: st_centroid assumes attributes are constant over geometries
#> Warning: neighbour object has 120 sub-graphs

boston_listw2 <- nb2listwdist(dnearneigh(st_centroid(boston.tr_utm), 5000, 9500), 
    boston.tr_utm, type = "dpd", alpha = 2, zero.policy = TRUE, dmax = 9500)
#> Warning: st_centroid assumes attributes are constant over geometries

boston_RM_gsi_1 <- localGS(boston.tr_utm, boston_listw1, 1, 2000, "RM", FALSE)
boston_RM_gsi_2 <- localGS(boston.tr_utm, boston_listw2, 2000, 9500, "RM", FALSE)
# }
```
