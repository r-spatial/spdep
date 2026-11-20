# Compute Normalized Moran's \\I\\

A simple function to compute the normalized "\\I_6\\" Moran's \\I\\
statistic proposed by Tillé et al. (2026), called by
`normalized.moran.test` and `normalized.moran.mc`.

## Usage

``` r
normalized.moran(
  x,
  listw,
  zero.policy = attr(listw, "zero.policy"),
  NAOK = FALSE
)
```

## Arguments

- x:

  a numeric vector the same length as the neighbours list in listw

- listw:

  a `listw` object of style 'W' created for example by `nb2listw`

- zero.policy:

  default `attr(listw, "zero.policy")` as set when `listw` was created,
  if attribute not set, use global option value; if TRUE assign zero to
  the lagged value of zones without neighbours, if FALSE assign NA

- NAOK:

  if 'TRUE' then any 'NA' or 'NaN' or 'Inf' values in x are passed on to
  the foreign function. If 'FALSE', the presence of 'NA' or 'NaN' or
  'Inf' values is regarded as an error.

## Value

The normalized "\\I_6\\" Moran's \\I\\ of the values in `x` is computed,
as a numeric vector of length one.

## Author

Diego Giuliani <diego.giuliani@unitn.it>

## References

Tillé Y., Giuliani D., Dickson M.M., Espa G. (2026) Revisiting the
Normalization of the Moran's \\I\\ index: A Correlation-based Approach
with Inference, *International Statistical Review*

## See also

[normalized.moran.test](https://r-spatial.github.io/spdep/reference/normalized.moran.test.md),
[normalized.moran.mc](https://r-spatial.github.io/spdep/reference/normalized.moran.mc.md)

## Examples

``` r
library(spdep)
data(oldcol)
col.W <- nb2listw(COL.nb, style="W")
crime <- COL.OLD$CRIME
normalized.moran(crime, col.W)
#> [1] 0.6908337
is.na(crime) <- sample(1:length(crime), 10)
normalized.moran(crime, col.W, NAOK=TRUE)
#> Warning: NAs in lagged values
#> [1] 0.449429
```
