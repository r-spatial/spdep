# Normalized Moran's \\I\\ test for spatial autocorrelation

Normalized "\\I_6\\" Moran's test for spatial autocorrelation proposed
by Tillé et al. (2026) using a spatial weights matrix in weights list
form. The test assumes asymptotic normality, and results may be checked
against those of `normalized.moran.mc` permutations.

## Usage

``` r
normalized.moran.test(
  x,
  listw,
  zero.policy = attr(listw, "zero.policy"),
  alternative = "greater",
  na.action = na.fail,
  spChk = NULL,
  adjust.n = TRUE
)
```

## Arguments

- x:

  a numeric vector the same length as the neighbours list in listw

- listw:

  a `listw` object created for example by `nb2listw`

- zero.policy:

  default `attr(listw, "zero.policy")` as set when `listw` was created,
  if attribute not set, use global option value; if TRUE assign zero to
  the lagged value of zones without neighbours, if FALSE assign NA

- alternative:

  a character string specifying the alternative hypothesis, must be one
  of `greater` (default), `less` or `two.sided`

- na.action:

  a function (default `na.fail`), can also be `na.omit` or
  `na.exclude` - in these cases the weights list will be subsetted to
  remove NAs in the data. It may be necessary to set `zero.policy` to
  TRUE because this subsetting may create no-neighbour observations.
  Note that only weights lists created without using the `glist`
  argument to `nb2listw` may be subsetted. If `na.pass` is used, zero is
  substituted for NA values in calculating the spatial lag

- spChk:

  should the data vector names be checked against the spatial objects
  for identity integrity, TRUE, or FALSE, default NULL to use
  [`get.spChkOption()`](https://r-spatial.github.io/spdep/reference/set.spChkOption.md)

- adjust.n:

  default TRUE, if FALSE the number of observations is not adjusted for
  no-neighbour observations, if TRUE, the number of observations is
  adjusted

## Value

A list with class `htest` containing the following components:

- `statistic`:

  the value of the standard deviate of normalized Moran's I.

- `p.value`:

  the p-value of the test.

- `estimate`:

  the value of the observed normalized Moran's I, its expectation and
  variance under normality.

- `alternative`:

  a character string describing the alternative hypothesis.

- `method`:

  a character string giving the assumption used for calculating the
  standard deviate.

- `data.name`:

  a character string giving the name(s) of the data.

## Author

Diego Giuliani <diego.giuliani@unitn.it>

## References

Tillé Y., Giuliani D., Dickson M.M., Espa G. (2026) Revisiting the
Normalization of the Moran's \\I\\ index: A Correlation-based Approach
with Inference, *International Statistical Review*

## See also

[normalized.moran](https://r-spatial.github.io/spdep/reference/normalized.moran.md),
[normalized.moran.mc](https://r-spatial.github.io/spdep/reference/normalized.moran.mc.md)

## Examples

``` r
library(spdep)
data(oldcol)
coords.OLD <- cbind(COL.OLD$X, COL.OLD$Y)
normalized.moran.test(COL.OLD$CRIME, nb2listw(COL.nb, style="W"))
#> 
#>  Normalized Moran I test under normality
#> 
#> data:  COL.OLD$CRIME 
#> weights: nb2listw(COL.nb, style = "W")    
#> 
#> Normalized Moran I statistic standard deviate = 3.6727, p-value =
#> 0.00012
#> alternative hypothesis: greater
#> sample estimates:
#> Normalized Moran I statistic                  Expectation 
#>                   0.69083370                  -0.04291698 
#>                     Variance 
#>                   0.03991333 
#> 
crime <- COL.OLD$CRIME
is.na(crime) <- sample(1:length(crime), 10)
res <- try(normalized.moran.test(crime, nb2listw(COL.nb, style="W"),
                                 na.action=na.fail))
#> Error in na.fail.default(x) : missing values in object
normalized.moran.test(crime, nb2listw(COL.nb, style="W"), zero.policy=TRUE,
                      na.action=na.omit)
#> Warning: subsetting caused increase in subgraph count
#> Error in solve.default(D): Lapack routine dgesv: system is exactly singular: U[7,7] = 0
normalized.moran.test(crime, nb2listw(COL.nb, style="W"), zero.policy=TRUE,
                      na.action=na.exclude)
#> Warning: subsetting caused increase in subgraph count
#> Error in solve.default(D): Lapack routine dgesv: system is exactly singular: U[7,7] = 0
normalized.moran.test(crime, nb2listw(COL.nb, style="W"), na.action=na.pass)
#> Warning: NAs in lagged values
#> 
#>  Normalized Moran I test under normality
#> 
#> data:  crime 
#> weights: nb2listw(COL.nb, style = "W")    
#> 
#> Normalized Moran I statistic standard deviate = 3.2963, p-value =
#> 0.0004898
#> alternative hypothesis: greater
#> sample estimates:
#> Normalized Moran I statistic                  Expectation 
#>                   0.61562834                  -0.04291698 
#>                     Variance 
#>                   0.03991333 
#> 
```
