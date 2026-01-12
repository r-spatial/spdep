# Calculate the local univariate join count

The univariate local join count statistic is used to identify clusters
of rarely occurring binary variables. The binary variable of interest
should occur less than half of the time.

## Usage

``` r
local_joincount_uni(
  fx,
  chosen,
  listw,
  alternative = "two.sided",
  nsim = 199,
  iseed = NULL,
  ties.method = "average",
  no_repeat_in_row=FALSE
)
```

## Arguments

- fx:

  a factor with two levels; use of an ordered factor is not well
  understood.

- chosen:

  a scalar character containing the level of `fx` that should be
  considered the observed value (1).

- listw:

  a listw object containing binary weights created, for example, with
  `nbwlistw(nb, style = "B")`

- alternative:

  default `"greater"`. One of `"less"` or `"greater"`.

- nsim:

  the number of conditional permutation simulations

- iseed:

  default NULL, used to set the seed; the output will only be
  reproducible if the count of CPU cores across which computation is
  distributed is the same

- ties.method:

  default `"average"` passed through to `rank`, can take values accepted
  by `rank`: `c("average", "first", "last", "random", "max", "min")`,
  see [`rank`](https://rdrr.io/r/base/rank.html)

- no_repeat_in_row:

  default `FALSE`, if `TRUE`, sample conditionally in each row without
  replacements to avoid duplicate values,
  <https://github.com/r-spatial/spdep/issues/124>

## Value

a `data.frame` with six columns `BB` (observed BB - neighbour same as
focus), `Pr()` (pseudo-p from punif rank), `sim_rank` (simulation rank
with current ties.method), `p_sim_pysal_ge` (pseudo-p fromPySAL esda,
greater than or equal to observed BB), `p_sim_pysal_gt` (pseudo-p
fromPySAL esda, greater than observed BB), `largereq` (count of
simulated values greater than or equal to observed BB) after folding,
`olarger` (count of simulated values greater than observed BB) before
folding, `olargereq` (count of simulated values greater than or equal to
observed BB) before folding,and number of rows equal to the length of
`x`.

## Details

The local join count statistic requires a binary weights list which can
be generated with `nb2listw(nb, style = "B")`. Additionally, ensure that
the binary variable of interest is rarely occurring in no more than half
of observations.

P-values are estimated using a conditional permutation approach. This
creates a reference distribution from which the observed statistic is
compared. For more see [Geoda
Glossary](https://raw.githubusercontent.com/GeoDaCenter/GeoDaCenter.github.io/refs/heads/update-1.22/glossary.html#ppvalue).

The pseudo-p-values returned by freestanding Geoda and PySAL esda
correspond to `"res_min$p_sim_pysal_ge"` and ranked pseudo-p-values with
`tied.method="min"` (equivalently `"last"`; `"res_min$p_sim_pysal_gt"`
corresponds to ranked pseudo-p-values with `tied.method="max"`
(equivalently `"first"`).

## References

Anselin, L., & Li, X. (2019). Operational Local Join Count Statistics
for Cluster Detection. Journal of geographical systems, 21(2), 189–210.
[doi:10.1007/s10109-019-00299-x](https://doi.org/10.1007/s10109-019-00299-x)

## Author

Josiah Parry <josiah.parry@gmail.com>

## Examples

``` r
data(oldcol)
fx <- as.factor(ifelse(COL.OLD$CRIME < 35, "low-crime", "high-crime"))
listw <- nb2listw(COL.nb, style = "B")
set.seed(1)
res_min <- local_joincount_uni(fx, chosen = "high-crime", listw, nsim=999,
 alternative="two.sided", ties.method="min")
cor(res_min[,2], res_min$p_sim_pysal_ge, use="complete.obs")
#> [1] 1
res_max <- local_joincount_uni(fx, chosen = "high-crime", listw, nsim=999,
 alternative="two.sided", ties.method="max")
cor(res_max[,2], res_max$p_sim_pysal_gt, use="complete.obs")
#> [1] 1
```
