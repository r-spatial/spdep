# Saddlepoint approximation of global Moran's I test

The function implements Tiefelsdorf's application of the Saddlepoint
approximation to global Moran's I's reference distribution.

## Usage

``` r
lm.morantest.sad(model, listw, zero.policy=attr(listw, "zero.policy"),
  alternative="greater", spChk=NULL, resfun=weighted.residuals,
  tol=.Machine$double.eps^0.5, maxiter=1000, tol.bounds=0.0001,
  zero.tol = 1e-07, Omega=NULL, save.M=NULL, save.U=NULL)
# S3 method for class 'moransad'
print(x, ...)
# S3 method for class 'moransad'
summary(object, ...)
# S3 method for class 'summary.moransad'
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

- tol:

  the desired accuracy (convergence tolerance) for `uniroot`

- maxiter:

  the maximum number of iterations for `uniroot`

- tol.bounds:

  offset from bounds for `uniroot`

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

- x:

  object to be printed

- object:

  object to be summarised

- ...:

  arguments to be passed through

## Details

The function involves finding the eigenvalues of an n by n matrix, and
numerically finding the root for the Saddlepoint approximation, and
should therefore only be used with care when n is large.

## Value

A list of class `moransad` with the following components:

- statistic:

  the value of the saddlepoint approximation of the standard deviate of
  global Moran's I.

- p.value:

  the p-value of the test.

- estimate:

  the value of the observed global Moran's I.

- alternative:

  a character string describing the alternative hypothesis.

- method:

  a character string giving the method used.

- data.name:

  a character string giving the name(s) of the data.

- internal1:

  Saddlepoint omega, r and u

- internal2:

  f.root, iter and estim.prec from `uniroot`

- df:

  degrees of freedom

- tau:

  eigenvalues (excluding zero values)

## References

Tiefelsdorf, M. 2002 The Saddlepoint approximation of Moran's I and
local Moran's Ii reference distributions and their numerical evaluation.
Geographical Analysis, 34, pp. 187–206; Bivand RS, Wong DWS 2018
Comparing implementations of global and local indicators of spatial
association. TEST, 27(3), 716–748
[doi:10.1007/s11749-018-0599-x](https://doi.org/10.1007/s11749-018-0599-x)

## Author

Roger Bivand <Roger.Bivand@nhh.no>

## See also

[`lm.morantest`](https://r-spatial.github.io/spdep/reference/lm.morantest.md)

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
summary(lm.morantest.sad(e.lm, nb2listw(eire.nb)))
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
#>  Expectation     Variance Std. deviate     Skewness     Kurtosis      Minimum 
#>  -0.05877741   0.01473183   3.25753938   0.31336881   3.05047361  -0.67545810 
#>      Maximum        omega        sad.r        sad.u 
#>   0.89091555   0.76549075   2.77616585   4.36854629 
#>       f.root         iter   estim.prec 
#> 6.996834e-14 1.100000e+01           NA 
e.wlm <- lm(OWNCONS ~ ROADACC, data=eire, weights=RETSALE)
lm.morantest(e.wlm, nb2listw(eire.nb), resfun=rstudent)
#> 
#>  Global Moran I for regression residuals
#> 
#> data:  
#> model: lm(formula = OWNCONS ~ ROADACC, data = eire, weights = RETSALE)
#> weights: nb2listw(eire.nb)
#> 
#> Moran I statistic standard deviate = 3.1385, p-value = 0.0008491
#> alternative hypothesis: greater
#> sample estimates:
#> Observed Moran I      Expectation         Variance 
#>       0.34500329      -0.04049313       0.01508687 
#> 
lm.morantest.sad(e.wlm, nb2listw(eire.nb), resfun=rstudent)
#> 
#>  Saddlepoint approximation for global Moran's I (Barndorff-Nielsen
#>  formula)
#> 
#> data:  
#> model:lm(formula = OWNCONS ~ ROADACC, data = eire, weights = RETSALE)
#> weights: nb2listw(eire.nb)
#> 
#> Saddlepoint approximation = 2.8708, p-value = 0.002047
#> alternative hypothesis: greater
#> sample estimates:
#> Observed Moran I 
#>        0.3450033 
#> 
```
