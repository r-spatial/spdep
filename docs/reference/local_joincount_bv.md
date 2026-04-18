# Calculate the local bivariate join count

The bivariate join count (BJC) evaluates event occurrences in predefined
regions and tests if the co-occurrence of events deviates from complete
spatial randomness.

## Usage

``` r
local_joincount_bv(
  x,
  z,
  listw,
  nsim = 199,
  alternative = "two.sided"
)
```

## Arguments

- x:

  a binary variable either numeric or logical

- z:

  a binary variable either numeric or logical with the same length as
  `x`

- listw:

  a listw object containing binary weights created, for example, with
  `nb2listw(nb, style = "B")`

- nsim:

  the number of conditional permutation simulations

- alternative:

  default `"greater"`. One of `"less"` or `"greater"`.

## Details

There are two cases that are evaluated in the bivariate join count. The
first being in-situ colocation (CLC) where xi = 1 and zi = 1. The second
is the general form of the bivariate join count (BJC) that is used when
there is no in-situ colocation.

The BJC case "is useful when x and z cannot occur in the same location,
such as when x and z correspond to two different values of a single
categorical variable" or "when x and z can co-locate, but do not"
(Anselin and Li, 2019). Whereas the CLC case is useful in evaluating
simultaneous occurrences of events.

The local bivariate join count statistic requires a binary weights list
which can be generated with `nb2listw(nb, style = "B")`.

P-values are only reported for those regions that match the CLC or BJC
criteria. Others will not have an associated p-value.

P-values are estimated using a conditional permutation approach. This
creates a reference distribution from which the observed statistic is
compared.

## Value

a `data.frame` with two columns `join_count` and `p_sim` and number of
rows equal to the length of arguments `x`.

## References

Anselin, L., & Li, X. (2019). Operational Local Join Count Statistics
for Cluster Detection. Journal of geographical systems, 21(2), 189–210.
[doi:10.1007/s10109-019-00299-x](https://doi.org/10.1007/s10109-019-00299-x)

## Author

Josiah Parry <josiah.parry@gmail.com>

## Examples

``` r
data("oldcol")
listw <- nb2listw(COL.nb, style = "B")
# Colocation case
x <- COL.OLD[["CP"]]
z <- COL.OLD[["EW"]]
set.seed(1)
res <- local_joincount_bv(x, z, listw)
na.omit(res)
#>    join_count p_sim
#> 9           2 0.480
#> 19          3 0.180
#> 20          4 0.020
#> 21          3 0.040
#> 22          3 0.060
#> 23          4 0.070
#> 24          5 0.005
#> 27          2 0.300
#> 28          3 0.050
#> 29          5 0.010
#> 30          5 0.010
#> 31          9 0.005
#> 32          4 0.090
#> 33          4 0.055
# no colocation case
z <- 1 - x
set.seed(1)
res <- local_joincount_bv(x, z, listw)
na.omit(res)
#>    join_count p_sim
#> 9           2 0.065
#> 19          3 0.300
#> 20          1 0.015
#> 21          1 0.040
#> 22          1 0.045
#> 23          3 0.165
#> 24          1 0.000
#> 27          2 0.230
#> 28          1 0.050
#> 29          1 0.000
#> 36          1 0.000
#> 37          2 0.005
#> 38          2 0.185
#> 39          2 0.265
#> 40          1 0.005
#> 43          1 0.070
```
