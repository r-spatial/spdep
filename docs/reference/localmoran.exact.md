# Exact local Moran's Ii tests

`localmoran.exact` provides exact local Moran's Ii tests under the null
hypothesis, while `localmoran.exact.alt` provides exact local Moran's Ii
tests under the alternative hypothesis. In this case, the model may be a
fitted model derived from a model fitted by
[`spatialreg::errorsarlm`](https://r-spatial.github.io/spatialreg/reference/ML_models.html),
with the covariance matrix can be passed through the `Omega=` argument.

## Usage

``` r
localmoran.exact(model, select, nb, glist = NULL, style = "W", 
 zero.policy = NULL, alternative = "two.sided", spChk = NULL, 
 resfun = weighted.residuals, save.Vi = FALSE, useTP=FALSE, truncErr=1e-6, 
 zeroTreat=0.1)
localmoran.exact.alt(model, select, nb, glist = NULL, style = "W",
 zero.policy = NULL, alternative = "two.sided", spChk = NULL,
 resfun = weighted.residuals, Omega = NULL, save.Vi = FALSE,
 save.M = FALSE, useTP=FALSE, truncErr=1e-6, zeroTreat=0.1)
# S3 method for class 'localmoranex'
print(x, ...)
# S3 method for class 'localmoranex'
as.data.frame(x, row.names=NULL, optional=FALSE, ...)
```

## Arguments

- model:

  an object of class `lm` returned by `lm` (assuming no global spatial
  autocorrelation), or an object of class `sarlm` returned by a spatial
  simultaneous autoregressive model fit (assuming global spatial
  autocorrelation represented by the model spatial coefficient); weights
  may be specified in the `lm` fit, but offsets should not be used

- select:

  an integer vector of the id. numbers of zones to be tested; if
  missing, all zones

- nb:

  a list of neighbours of class `nb`

- glist:

  a list of general weights corresponding to neighbours

- style:

  can take values W, B, C, and S

- zero.policy:

  default NULL, use global option value; if TRUE assign zero to the
  lagged value of zones without neighbours, if FALSE assign NA

- alternative:

  a character string specifying the alternative hypothesis, must be one
  of greater (default), less or two.sided.

- spChk:

  should the data vector names be checked against the spatial objects
  for identity integrity, TRUE, or FALSE, default NULL to use
  [`get.spChkOption()`](https://r-spatial.github.io/spdep/reference/set.spChkOption.md)

- resfun:

  default: weighted.residuals; the function to be used to extract
  residuals from the `lm` object, may be `residuals`,
  `weighted.residuals`, `rstandard`, or `rstudent`

- Omega:

  A SAR process matrix may be passed in to test an alternative
  hypothesis, for example
  `Omega <- invIrW(listw, rho=0.1); Omega <- tcrossprod(Omega)`,
  [`chol()`](https://rdrr.io/pkg/Matrix/man/chol-methods.html) is taken
  internally

- save.Vi:

  if TRUE, return the star-shaped weights lists for each zone tested

- save.M:

  if TRUE, save a list of left and right M products

- useTP:

  default FALSE, if TRUE, use truncation point in integration rather
  than upper=Inf, see Tiefelsdorf (2000), eq. 6.7, p.69

- truncErr:

  when useTP=TRUE, pass truncation error to truncation point function

- zeroTreat:

  when useTP=TRUE, pass zero adjustment to truncation point function

- x:

  object to be printed

- row.names:

  ignored argument to `as.data.frame.localmoranex`; row names assigned
  from localmoranex object

- optional:

  ignored argument to `as.data.frame.localmoranex`; row names assigned
  from localmoranex object

- ...:

  arguments to be passed through

## Value

A list with class `localmoranex` containing "select" lists, each with
class `moranex` with the following components:

- statistic:

  the value of the exact standard deviate of global Moran's I.

- p.value:

  the p-value of the test.

- estimate:

  the value of the observed local Moran's Ii.

- method:

  a character string giving the method used.

- alternative:

  a character string describing the alternative hypothesis.

- gamma:

  eigenvalues (two extreme values for null, vector for alternative)

- oType:

  usually set to "E", but set to "N" if the integration leads to an out
  of domain value for `qnorm`, when the Normal assumption is
  substituted. This only occurs when the output p-value would be very
  close to zero

- data.name:

  a character string giving the name(s) of the data.

- df:

  degrees of freedom

- i:

  zone tested

- Vi:

  zone tested

When the alternative is being tested, a list of left and right M
products in attribute M.

## References

Bivand RS, Müller W, Reder M (2009) Power calculations for global and
local Moran’s I. Comput Stat Data Anal 53:2859–2872; Bivand RS, Wong DWS
2018 Comparing implementations of global and local indicators of spatial
association. TEST, 27(3), 716–748
[doi:10.1007/s11749-018-0599-x](https://doi.org/10.1007/s11749-018-0599-x)

## Author

Markus Reder and Roger Bivand

## See also

[`lm.morantest.exact`](https://r-spatial.github.io/spdep/reference/lm.morantest.exact.md),
[`localmoran.sad`](https://r-spatial.github.io/spdep/reference/localmoran.sad.md)

## Examples

``` r
eire <- st_read(system.file("shapes/eire.gpkg", package="spData")[1])
#> Reading layer `eire' from data source 
#>   `/home/rsb/lib/r_libs/spData/shapes/eire.gpkg' using driver `GPKG'
#> Simple feature collection with 26 features and 10 fields
#> Geometry type: MULTIPOLYGON
#> Dimension:     XY
#> Bounding box:  xmin: -4.12 ymin: 5768 xmax: 300.82 ymax: 6119.25
#> Projected CRS: Undefined Cartesian SRS with unknown unit
row.names(eire) <- as.character(eire$names)
eire.nb <- poly2nb(eire)
e.lm <- lm(OWNCONS ~ ROADACC, data=eire)
localmoran.sad(e.lm, nb=eire.nb)
#>              Local Morans I Saddlepoint    Pr. (Sad)
#> 1 Carlow         0.21699668  0.95074844 3.417321e-01
#> 2 Cavan         -0.37257361 -1.00603119 3.144006e-01
#> 3 Clare          0.23197510  0.67166518 5.017969e-01
#> 4 Cork           0.78193548  1.74761575 8.053059e-02
#> 5 Donegal       -1.69064059 -1.72031078 8.537596e-02
#> 6 Dublin        -0.16069692 -0.35212627 7.247436e-01
#> 7 Galway         1.31371473  2.66849536 7.619183e-03
#> 8 Kerry          0.36534866  0.78073279 4.349597e-01
#> 9 Kildare       -0.02557544  0.04167665 9.667565e-01
#> 10 Kilkenny      0.57684331  1.70897697 8.745521e-02
#> 11 Laoghis      -0.05951798 -0.12155465 9.032517e-01
#> 12 Leitrim       0.38484587  1.47227033 1.409479e-01
#> 13 Limerick      0.11817987  0.45727712 6.474719e-01
#> 14 Longford      1.41643200  2.51113769 1.203427e-02
#> 15 Louth         0.56242920  1.07441571 2.826364e-01
#> 16 Mayo          0.87572704  2.05251226 4.011990e-02
#> 17 Meath         0.00367856  0.12813539 8.980418e-01
#> 18 Monaghan      0.55098311  1.23999193 2.149784e-01
#> 19 Offaly        0.15155556  0.80786519 4.191682e-01
#> 20 Roscommon     2.04368839  4.53187292 5.846302e-06
#> 21 Sligo        -0.47579871 -0.94578114 3.442602e-01
#> 22 Tipperary    -0.03454106 -0.06919691 9.448329e-01
#> 23 Waterford     0.85723423  1.91385108 5.563919e-02
#> 24 Westmeath     0.45138572  1.36017204 1.737755e-01
#> 25 Wexford       0.64371834  1.63188492 1.027037e-01
#> 26 Wicklow       0.02441950  0.21197000 8.321304e-01
localmoran.exact(e.lm, nb=eire.nb)
#>              Local Morans I    Exact SD  Pr. (exact)
#> 1 Carlow         0.21699668  1.02706083 3.043918e-01
#> 2 Cavan         -0.37257361 -1.04864802 2.943401e-01
#> 3 Clare          0.23197510  0.76362894 4.450884e-01
#> 4 Cork           0.78193548  1.77727656 7.552275e-02
#> 5 Donegal       -1.69064059 -1.74428756 8.110896e-02
#> 6 Dublin        -0.16069692 -0.44236119 6.582279e-01
#> 7 Galway         1.31371473  2.69199974 7.102500e-03
#> 8 Kerry          0.36534866  0.85742696 3.912090e-01
#> 9 Kildare       -0.02557544 -0.03475476 9.722753e-01
#> 10 Kilkenny      0.57684331  1.74146177 8.160267e-02
#> 11 Laoghis      -0.05951798 -0.21824035 8.272419e-01
#> 12 Leitrim       0.38484587  1.51434641 1.299380e-01
#> 13 Limerick      0.11817987  0.56922630 5.692026e-01
#> 14 Longford      1.41643200  2.53491837 1.124735e-02
#> 15 Louth         0.56242920  1.12775107 2.594251e-01
#> 16 Mayo          0.87572704  2.08125803 3.741029e-02
#> 17 Meath         0.00367856  0.16372685 8.699462e-01
#> 18 Monaghan      0.55098311  1.28435459 1.990179e-01
#> 19 Offaly        0.15155556  0.89537870 3.705847e-01
#> 20 Roscommon     2.04368839  4.55244870 5.302509e-06
#> 21 Sligo        -0.47579871 -0.98101752 3.265841e-01
#> 22 Tipperary    -0.03454106 -0.16132608 8.718366e-01
#> 23 Waterford     0.85723423  1.94188723 5.215075e-02
#> 24 Westmeath     0.45138572  1.40091422 1.612397e-01
#> 25 Wexford       0.64371834  1.66488051 9.593660e-02
#> 26 Wicklow       0.02441950  0.29717701 7.663314e-01
localmoran.exact(e.lm, nb=eire.nb, useTP=TRUE)
#>              Local Morans I    Exact SD  Pr. (exact)
#> 1 Carlow         0.21699668  1.02706127 3.043916e-01
#> 2 Cavan         -0.37257361 -1.04864834 2.943400e-01
#> 3 Clare          0.23197510  0.76362936 4.450881e-01
#> 4 Cork           0.78193548  1.77727666 7.552273e-02
#> 5 Donegal       -1.69064059 -1.74428741 8.110899e-02
#> 6 Dublin        -0.16069692 -0.44236158 6.582276e-01
#> 7 Galway         1.31371473  2.69199162 7.102673e-03
#> 8 Kerry          0.36534866  0.85742736 3.912087e-01
#> 9 Kildare       -0.02557544 -0.03540278 9.717586e-01
#> 10 Kilkenny      0.57684331  1.74146192 8.160264e-02
#> 11 Laoghis      -0.05951798 -0.21823772 8.272439e-01
#> 12 Leitrim       0.38484587  1.51434680 1.299379e-01
#> 13 Limerick      0.11817987  0.56922674 5.692023e-01
#> 14 Longford      1.41643200  2.53491383 1.124750e-02
#> 15 Louth         0.56242920  1.12775145 2.594249e-01
#> 16 Mayo          0.87572704  2.08125778 3.741032e-02
#> 17 Meath         0.00367856  0.15714487 8.751307e-01
#> 18 Monaghan      0.55098311  1.28435495 1.990178e-01
#> 19 Offaly        0.15155556  0.89537914 3.705844e-01
#> 20 Roscommon     2.04368839  4.53362873 5.797890e-06
#> 21 Sligo        -0.47579871 -0.98101780 3.265840e-01
#> 22 Tipperary    -0.03454106 -0.16116590 8.719627e-01
#> 23 Waterford     0.85723423  1.94188711 5.215077e-02
#> 24 Westmeath     0.45138572  1.40091457 1.612396e-01
#> 25 Wexford       0.64371834  1.66488069 9.593656e-02
#> 26 Wicklow       0.02441950  0.29790569 7.657751e-01
run <- FALSE
if (requireNamespace("spatialreg", quietly=TRUE)) run <- TRUE
if (run) {
e.errorsar <- spatialreg::errorsarlm(OWNCONS ~ ROADACC, data=eire,
 listw=nb2listw(eire.nb))
lm.target <- lm(e.errorsar$tary ~ e.errorsar$tarX - 1)
localmoran.exact.alt(lm.target, nb=eire.nb)
}
#>              Local Morans I    Exact SD Pr. (exact)
#> 1 Carlow         0.17958462  0.83313548 0.404768327
#> 2 Cavan         -0.24752628 -0.82861508 0.407322251
#> 3 Clare         -0.27901334 -0.72313220 0.469598670
#> 4 Cork           0.37808655  1.15083938 0.249798320
#> 5 Donegal       -1.01894688 -1.50937607 0.131202704
#> 6 Dublin        -0.18171297 -0.50995181 0.610085224
#> 7 Galway         1.02193390  2.27605549 0.022842689
#> 8 Kerry         -0.94967914 -1.46954733 0.141684395
#> 9 Kildare        0.07005053  0.56174770 0.574287931
#> 10 Kilkenny      0.43022231  1.39078833 0.164289622
#> 11 Laoghis      -0.12239133 -0.45750046 0.647311369
#> 12 Leitrim      -0.24203970 -0.84925558 0.395739092
#> 13 Limerick     -0.03214546 -0.13098232 0.895789292
#> 14 Longford      0.38307454  1.17909760 0.238359315
#> 15 Louth         0.21301552  0.62068302 0.534808204
#> 16 Mayo          0.93971200  1.94079901 0.052282661
#> 17 Meath         0.12484415  0.78747315 0.431004937
#> 18 Monaghan     -0.16109919 -0.48277561 0.629255068
#> 19 Offaly       -0.00632492  0.03008421 0.975999892
#> 20 Roscommon     1.02089429  2.64333471 0.008209384
#> 21 Sligo        -2.01629233 -2.54211789 0.011018300
#> 22 Tipperary    -0.10810709 -0.51063241 0.609608478
#> 23 Waterford     0.44099279  1.28152280 0.200010096
#> 24 Westmeath    -0.06329661 -0.26575204 0.790430182
#> 25 Wexford       0.30764883  1.02115599 0.307180540
#> 26 Wicklow      -0.01696406 -0.04109885 0.967217095
if (run) {
Omega <- spatialreg::invIrW(nb2listw(eire.nb), rho=e.errorsar$lambda)
Omega1 <- tcrossprod(Omega)
localmoran.exact.alt(lm.target, nb=eire.nb, Omega=Omega1)
}
#>              Local Morans I     Exact SD  Pr. (exact)
#> 1 Carlow         0.17958462  0.009413917 0.9924888917
#> 2 Cavan         -0.24752628 -1.416786314 0.1565454123
#> 3 Clare         -0.27901334 -1.512813055 0.1303271518
#> 4 Cork           0.37808655  0.078686494 0.9372819882
#> 5 Donegal       -1.01894688 -1.935314161 0.0529517556
#> 6 Dublin        -0.18171297 -1.226543250 0.2199943101
#> 7 Galway         1.02193390  1.477322295 0.1395892645
#> 8 Kerry         -0.94967914 -2.427948363 0.0151845050
#> 9 Kildare        0.07005053 -0.042522018 0.9660825597
#> 10 Kilkenny      0.43022231  0.488471215 0.6252161095
#> 11 Laoghis      -0.12239133 -1.091910862 0.2748722828
#> 12 Leitrim      -0.24203970 -1.631825829 0.1027161825
#> 13 Limerick     -0.03214546 -0.921751556 0.3566581814
#> 14 Longford      0.38307454  0.183273193 0.8545836668
#> 15 Louth         0.21301552 -0.184065307 0.8539622228
#> 16 Mayo          0.93971200  1.182398032 0.2370478039
#> 17 Meath         0.12484415 -0.060499926 0.9517574751
#> 18 Monaghan     -0.16109919 -1.342787508 0.1793407894
#> 19 Offaly       -0.00632492 -0.409598788 0.6821002864
#> 20 Roscommon     1.02089429  1.238387417 0.2155724439
#> 21 Sligo        -2.01629233 -3.623597237 0.0002905339
#> 22 Tipperary    -0.10810709 -1.227333133 0.2196974073
#> 23 Waterford     0.44099279  0.433358039 0.6647546795
#> 24 Westmeath    -0.06329661 -0.956519228 0.3388099749
#> 25 Wexford       0.30764883  0.141381249 0.8875687669
#> 26 Wicklow      -0.01696406 -0.783977878 0.4330530912
if (run) {
localmoran.exact.alt(lm.target, nb=eire.nb, Omega=Omega1, useTP=TRUE)
}
#>              Local Morans I     Exact SD  Pr. (exact)
#> 1 Carlow         0.17958462  0.009414171 0.9924886890
#> 2 Cavan         -0.24752628 -1.416786809 0.1565452674
#> 3 Clare         -0.27901334 -1.512813519 0.1303270339
#> 4 Cork           0.37808655  0.078686628 0.9372818820
#> 5 Donegal       -1.01894688 -1.935314436 0.0529517219
#> 6 Dublin        -0.18171297 -1.226543734 0.2199941280
#> 7 Galway         1.02193390  1.477322399 0.1395892366
#> 8 Kerry         -0.94967914 -2.427949220 0.0151844692
#> 9 Kildare        0.07005053 -0.042526750 0.9660787872
#> 10 Kilkenny      0.43022231  0.488471349 0.6252160151
#> 11 Laoghis      -0.12239133 -1.091912026 0.2748717712
#> 12 Leitrim      -0.24203970 -1.631826543 0.1027160320
#> 13 Limerick     -0.03214546 -0.921785589 0.3566404260
#> 14 Longford      0.38307454  0.183273325 0.8545835633
#> 15 Louth         0.21301552 -0.184065128 0.8539623632
#> 16 Mayo          0.93971200  1.182398127 0.2370477660
#> 17 Meath         0.12484415 -0.060499243 0.9517580191
#> 18 Monaghan     -0.16109919 -1.342788191 0.1793405683
#> 19 Offaly       -0.00632492 -0.407816672 0.6834082724
#> 20 Roscommon     1.02089429  1.238387488 0.2155724176
#> 21 Sligo        -2.01629233 -3.623590126 0.0002905419
#> 22 Tipperary    -0.10810709 -1.227335007 0.2196967034
#> 23 Waterford     0.44099279  0.433358168 0.6647545856
#> 24 Westmeath    -0.06329661 -0.956505669 0.3388168216
#> 25 Wexford       0.30764883  0.141381402 0.8875686458
#> 26 Wicklow      -0.01696406 -0.784451760 0.4327750770
```
