# Permutation test for same colour join count statistics

A permutation test for same colour join count statistics calculated by
using nsim random permutations of fx for the given spatial weighting
scheme, to establish the ranks of the observed statistics (for each
colour) in relation to the nsim simulated values.

## Usage

``` r
joincount.mc(fx, listw, nsim, zero.policy=attr(listw, "zero.policy"),
 alternative="greater", spChk=NULL)
```

## Arguments

- fx:

  a factor of the same length as the neighbours and weights objects in
  listw; use of an ordered factor is not well understood

- listw:

  a `listw` object created for example by `nb2listw`

- nsim:

  number of permutations

- zero.policy:

  default `attr(listw, "zero.policy")` as set when `listw` was created,
  if attribute not set, use global option value; if TRUE assign zero to
  the lagged value of zones without neighbours, if FALSE assign NA

- alternative:

  a character string specifying the alternative hypothesis, must be one
  of "greater" (default), "two.sided", or "less".

- spChk:

  should the data vector names be checked against the spatial objects
  for identity integrity, TRUE, or FALSE, default NULL to use
  [`get.spChkOption()`](https://r-spatial.github.io/spdep/reference/set.spChkOption.md)

## Value

A list with class `jclist` of lists with class `htest` and `mc.sim` for
each of the k colours containing the following components:

- statistic:

  the value of the observed statistic.

- parameter:

  the rank of the observed statistic.

- method:

  a character string giving the method used.

- data.name:

  a character string giving the name(s) of the data.

- p.value:

  the pseudo p-value of the test.

- alternative:

  a character string describing the alternative hypothesis.

- estimate:

  the mean and variance of the simulated distribution.

- res:

  nsim simulated values of statistic, the final element is the observed
  statistic

## References

Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 63-5.

## Author

Roger Bivand <Roger.Bivand@nhh.no>

## See also

[`joincount.test`](https://r-spatial.github.io/spdep/reference/joincount.test.md)

## Examples

``` r
data(oldcol)
HICRIME <- cut(COL.OLD$CRIME, breaks=c(0,35,80), labels=c("low","high"))
names(HICRIME) <- rownames(COL.OLD)
joincount.mc(HICRIME, nb2listw(COL.nb, style="B"), nsim=99, alternative="two.sided")
#> 
#>  Monte-Carlo simulation of join-count statistic
#> 
#> data:  HICRIME 
#> weights: nb2listw(COL.nb, style = "B") 
#> number of simulations + 1: 100 
#> 
#> Join-count statistic for low = 34, rank of observed statistic = 89.5,
#> p-value = 0.21
#> alternative hypothesis: two.sided
#> sample estimates:
#>     mean of simulation variance of simulation 
#>               28.92929               17.78066 
#> 
#> 
#>  Monte-Carlo simulation of join-count statistic
#> 
#> data:  HICRIME 
#> weights: nb2listw(COL.nb, style = "B") 
#> number of simulations + 1: 100 
#> 
#> Join-count statistic for high = 54, rank of observed statistic = 100,
#> p-value < 2.2e-16
#> alternative hypothesis: two.sided
#> sample estimates:
#>     mean of simulation variance of simulation 
#>               26.86869               20.11523 
#> 
joincount.test(HICRIME, nb2listw(COL.nb, style="B"), alternative="two.sided")
#> 
#>  Join count test under nonfree sampling
#> 
#> data:  HICRIME 
#> weights: nb2listw(COL.nb, style = "B") 
#> 
#> Std. deviate for low = 1.0141, p-value = 0.3105
#> alternative hypothesis: two.sided
#> sample estimates:
#> Same colour statistic           Expectation              Variance 
#>              34.00000              29.59184              18.89550 
#> 
#> 
#>  Join count test under nonfree sampling
#> 
#> data:  HICRIME 
#> weights: nb2listw(COL.nb, style = "B") 
#> 
#> Std. deviate for high = 6.3307, p-value = 2.44e-10
#> alternative hypothesis: two.sided
#> sample estimates:
#> Same colour statistic           Expectation              Variance 
#>              54.00000              27.22449              17.88838 
#> 
```
