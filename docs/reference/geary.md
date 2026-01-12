# Compute Geary's C

A simple function to compute Geary's C, called by `geary.test` and
`geary.mc`; \$\$C = \frac{(n-1)}{2\sum\_{i=1}^{n}\sum\_{j=1}^{n}w\_{ij}}
\frac{\sum\_{i=1}^{n}\sum\_{j=1}^{n}w\_{ij}(x_i-x_j)^2}{\sum\_{i=1}^{n}(x_i -
\bar{x})^2} \$\$ `geary.intern` is an internal function used to vary the
similarity criterion.

## Usage

    geary(x, listw, n, n1, S0, zero.policy=attr(listw, "zero.policy"), scale=TRUE)
    <!-- %geary.intern(x, listw, n, zero.policy, type="geary") -->

## Arguments

- x:

  a numeric vector the same length as the neighbours list in listw

- listw:

  a `listw` object created for example by `nb2listw`

- n:

  number of zones

- n1:

  n - 1

- S0:

  global sum of weights

- zero.policy:

  default `attr(listw, "zero.policy")` as set when `listw` was created,
  if attribute not set, use global option value; if TRUE assign zero to
  the lagged value of zones without neighbours, if FALSE assign NA

&nbsp;

- scale:

  default TRUE, may be FALSE to revert changes made to accommodate
  `localC` in November 2021 (see \#151)

## Value

a list with

- C:

  Geary's C

- K:

  sample kurtosis of x

## References

Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 17.

## Author

Roger Bivand <Roger.Bivand@nhh.no>

## See also

[`geary.test`](https://r-spatial.github.io/spdep/reference/geary.test.md),
[`geary.mc`](https://r-spatial.github.io/spdep/reference/geary.mc.md),
[`sp.mantel.mc`](https://r-spatial.github.io/spdep/reference/sp.mantel.mc.md)

## Examples

``` r
data(oldcol)
col.W <- nb2listw(COL.nb, style="W")
str(geary(COL.OLD$CRIME, col.W, length(COL.nb), length(COL.nb)-1,
 Szero(col.W)))
#> List of 2
#>  $ C: num 0.53
#>  $ K: num 2.23
```
