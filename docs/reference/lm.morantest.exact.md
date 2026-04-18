# Exact global Moran's I test

The function implements Tiefelsdorf's exact global Moran's I test.

## Usage

``` r
lm.morantest.exact(model, listw, zero.policy = attr(listw, "zero.policy"),
 alternative = "greater", spChk = NULL, resfun = weighted.residuals,
 zero.tol = 1e-07, Omega=NULL, save.M=NULL, save.U=NULL, useTP=FALSE,
 truncErr=1e-6, zeroTreat=0.1)
# S3 method for class 'moranex'
print(x, ...)
```

## Arguments

- model:

  an object of class `lm` returned by `lm`; weights may be specified in
  the `lm` fit, but offsets should not be used

- listw:

  a `listw` object created for example by `nb2listw`

- zero.policy:

  default `attr(listw, "zero.policy")` as set when `listw` was created,
  if attribute not set, use global option value; if TRUE assign zero to
  the lagged value of zones without neighbours, if FALSE assign NA

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

- zero.tol:

  tolerance used to find eigenvalues close to absolute zero

- Omega:

  A SAR process matrix may be passed in to test an alternative
  hypothesis, for example
  `Omega <- invIrW(listw, rho=0.1); Omega <- tcrossprod(Omega)`,
  [`chol()`](https://rdrr.io/pkg/Matrix/man/chol-methods.html) is taken
  internally

- save.M:

  return the full M matrix for use in `spdep:::exactMoranAlt`

- save.U:

  return the full U matrix for use in `spdep:::exactMoranAlt`

- useTP:

  default FALSE, if TRUE, use truncation point in integration rather
  than upper=Inf, see Tiefelsdorf (2000), eq. 6.7, p.69

- truncErr:

  when useTP=TRUE, pass truncation error to truncation point function

- zeroTreat:

  when useTP=TRUE, pass zero adjustment to truncation point function

- x:

  a moranex object

- ...:

  arguments to be passed through

## Value

A list of class `moranex` with the following components:

- statistic:

  the value of the exact standard deviate of global Moran's I.

- p.value:

  the p-value of the test.

- estimate:

  the value of the observed global Moran's I.

- method:

  a character string giving the method used.

- alternative:

  a character string describing the alternative hypothesis.

- gamma:

  eigenvalues (excluding zero values)

- oType:

  usually set to "E"

- data.name:

  a character string giving the name(s) of the data.

- df:

  degrees of freedom

## Author

Markus Reder and Roger Bivand

## References

Roger Bivand, Werner G. Müller and Markus Reder (2009) "Power
calculations for global and local Moran's I." *Computational Statistics
& Data Analysis* 53, 2859-2872.

## See also

[`lm.morantest.sad`](https://r-spatial.github.io/spdep/reference/lm.morantest.sad.md)

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
lm.morantest(e.lm, nb2listw(eire.nb))
#> 
#>  Global Moran I for regression residuals
#> 
#> data:  
#> model: lm(formula = OWNCONS ~ ROADACC, data = eire)
#> weights: nb2listw(eire.nb)
#> 
#> Moran I statistic standard deviate = 3.2575, p-value = 0.0005619
#> alternative hypothesis: greater
#> sample estimates:
#> Observed Moran I      Expectation         Variance 
#>       0.33660565      -0.05877741       0.01473183 
#> 
lm.morantest.sad(e.lm, nb2listw(eire.nb))
#> 
#>  Saddlepoint approximation for global Moran's I (Barndorff-Nielsen
#>  formula)
#> 
#> data:  
#> model:lm(formula = OWNCONS ~ ROADACC, data = eire)
#> weights: nb2listw(eire.nb)
#> 
#> Saddlepoint approximation = 2.9395, p-value = 0.001644
#> alternative hypothesis: greater
#> sample estimates:
#> Observed Moran I 
#>        0.3366057 
#> 
lm.morantest.exact(e.lm, nb2listw(eire.nb))
#> 
#>  Global Moran I statistic with exact p-value
#> 
#> data:  
#> model:lm(formula = OWNCONS ~ ROADACC, data = eire)
#> weights: nb2listw(eire.nb)
#> 
#> Exact standard deviate = 2.9316, p-value = 0.001686
#> alternative hypothesis: greater
#> sample estimates:
#> [1] 0.3366057
#> 
lm.morantest.exact(e.lm, nb2listw(eire.nb), useTP=TRUE)
#> 
#>  Global Moran I statistic with exact p-value
#> 
#> data:  
#> model:lm(formula = OWNCONS ~ ROADACC, data = eire)
#> weights: nb2listw(eire.nb)
#> 
#> Exact standard deviate = 2.9315, p-value = 0.001686
#> alternative hypothesis: greater
#> sample estimates:
#> [1] 0.3366057
#> 
```
