# Compute Moran's I

A simple function to compute Moran's I, called by `moran.test` and
`moran.mc`; \$\$I = \frac{n}{\sum\_{i=1}^{n}\sum\_{j=1}^{n}w\_{ij}}
\frac{\sum\_{i=1}^{n}\sum\_{j=1}^{n}w\_{ij}(x_i-\bar{x})(x_j-\bar{x})}{\sum\_{i=1}^{n}(x_i -
\bar{x})^2} \$\$

## Usage

``` r
moran(x, listw, n, S0, zero.policy=attr(listw, "zero.policy"), NAOK=FALSE)
```

## Arguments

- x:

  a numeric vector the same length as the neighbours list in listw

- listw:

  a `listw` object created for example by `nb2listw`

- n:

  number of zones

- S0:

  global sum of weights

- zero.policy:

  default `attr(listw, "zero.policy")` as set when `listw` was created,
  if attribute not set, use global option value; if TRUE assign zero to
  the lagged value of zones without neighbours, if FALSE assign NA

- NAOK:

  if 'TRUE' then any 'NA' or 'NaN' or 'Inf' values in x are passed on to
  the foreign function. If 'FALSE', the presence of 'NA' or 'NaN' or
  'Inf' values is regarded as an error.

## Value

a list of

- I:

  Moran's I

- K:

  sample kurtosis of x

## References

Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 17.

## Author

Roger Bivand <Roger.Bivand@nhh.no>

## See also

[`moran.test`](https://r-spatial.github.io/spdep/reference/moran.test.md),
[`moran.mc`](https://r-spatial.github.io/spdep/reference/moran.mc.md)

## Examples

``` r
data(oldcol)
col.W <- nb2listw(COL.nb, style="W")
crime <- COL.OLD$CRIME
str(moran(crime, col.W, length(COL.nb), Szero(col.W)))
#> List of 2
#>  $ I: num 0.511
#>  $ K: num 2.23
is.na(crime) <- sample(1:length(crime), 10)
str(moran(crime, col.W, length(COL.nb), Szero(col.W), NAOK=TRUE))
#> Warning: NAs in lagged values
#> List of 2
#>  $ I: num 0.283
#>  $ K: num 2.67
```
