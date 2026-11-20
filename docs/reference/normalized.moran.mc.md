# Permutation test for normalized Moran's \\I\\ statistic

A permutation test for the "\\I_6\\" normalized Moran's \\I\\ statistic
by Tillé et al. (2026) calculated by using `nsim` random permutations of
`x` for the given spatial weighting scheme, to establish the rank of the
observed statistic in relation to the `nsim` simulated values.

## Usage

``` r
normalized.moran.mc(
  x,
  listw,
  nsim,
  zero.policy = attr(listw, "zero.policy"),
  alternative = "greater",
  na.action = na.fail,
  spChk = NULL,
  return_boot = FALSE,
  adjust.n = TRUE
)
```

## Arguments

- x:

  a numeric vector the same length as the neighbours list in `listw`

- listw:

  a `listw` object created for example by nb2listw

- nsim:

  number of permutations

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

- return_boot:

  return an object of class `boot` from the equivalent permutation
  bootstrap rather than an object of class `htest`

- adjust.n:

  default TRUE, if FALSE the number of observations is not adjusted for
  no-neighbour observations, if TRUE, the number of observations is
  adjusted

## Value

A list with class `htest` and `mc.sim` containing the following
components:

- `statistic`:

  the value of the observed normalized Moran's I.

- `parameter`:

  the rank of the observed normalized Moran's I.

- `p.value`:

  the pseudo p-value of the test.

- `alternative`:

  a character string describing the alternative hypothesis.

- `method`:

  a character string giving the method used.

- `data.name`:

  a character string giving the name(s) of the data.

- `res`:

  `nsim` simulated values of statistic, final value is observed
  statistic.

## Author

Diego Giuliani <diego.giuliani@unitn.it>

## References

Tillé Y., Giuliani D., Dickson M.M., Espa G. (2026) Revisiting the
Normalization of the Moran's \\I\\ index: A Correlation-based Approach
with Inference, *International Statistical Review*

## See also

[normalized.moran](https://r-spatial.github.io/spdep/reference/normalized.moran.md),
[normalized.moran.test](https://r-spatial.github.io/spdep/reference/normalized.moran.test.md)

## Examples

``` r
library(spdep)
data(oldcol)
colw <- nb2listw(COL.nb, style="W")
nsim <- 99
set.seed(1234)
sim1 <- normalized.moran.mc(COL.OLD$CRIME, listw=colw, nsim=nsim)
sim1
#> 
#>  Monte-Carlo simulation of Normalized Moran I
#> 
#> data:  COL.OLD$CRIME 
#> weights: colw  
#> number of simulations + 1: 100 
#> 
#> statistic = 0.69083, observed rank = 100, p-value = 0.01
#> alternative hypothesis: greater
#> 
mean(sim1$res[1:nsim])
#> [1] -0.01658006
var(sim1$res[1:nsim])
#> [1] 0.03477783
summary(sim1$res[1:nsim])
#>     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#> -0.47987 -0.16358 -0.02060 -0.01658  0.12946  0.37872 
colold.lags <- nblag(COL.nb, 3)
set.seed(1234)
sim2 <- normalized.moran.mc(COL.OLD$CRIME, nb2listw(colold.lags[[2]],
                                                    style="W"), nsim=nsim)
summary(sim2$res[1:nsim])
#>     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#> -0.42877 -0.20490 -0.04378 -0.05787  0.08707  0.34391 
sim3 <- normalized.moran.mc(COL.OLD$CRIME, nb2listw(colold.lags[[3]],
                                                    style="W"), nsim=nsim)
summary(sim3$res[1:nsim])
#>     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#> -0.47870 -0.16714 -0.06838 -0.07033  0.01565  0.32550 
crime <- COL.OLD$CRIME
is.na(crime) <- sample(1:length(crime), 10)
try(normalized.moran.mc(crime, nb2listw(COL.nb, style="W"), nsim=99,
                        na.action=na.fail))
#> Error in na.fail.default(x) : missing values in object
normalized.moran.mc(crime, nb2listw(COL.nb, style="W"), nsim=99, zero.policy=TRUE,
                    na.action=na.omit)
#> 
#>  Monte-Carlo simulation of Normalized Moran I
#> 
#> data:  crime 
#> weights: nb2listw(COL.nb, style = "W") 
#> omitted: 5, 7, 9, 15, 19, 33, 40, 45, 47, 48 
#> number of simulations + 1: 100 
#> 
#> statistic = 0.60616, observed rank = 100, p-value = 0.01
#> alternative hypothesis: greater
#> 
normalized.moran.mc(crime, nb2listw(COL.nb, style="W"), nsim=99, zero.policy=TRUE,
                    return_boot=TRUE, na.action=na.omit)
#> NA observations omitted: 5, 7, 9, 15, 19, 33, 40, 45, 47, 48
#> 
#> DATA PERMUTATION
#> 
#> 
#> Call:
#> boot::boot(data = x, statistic = moran_boot, R = nsim, sim = "permutation", 
#>     listw = listw, zero.policy = zero.policy, parallel = parallel, 
#>     ncpus = ncpus, cl = cl)
#> 
#> 
#> Bootstrap Statistics :
#>      original     bias    std. error
#> t1* 0.6061552 -0.6615062   0.2414765
normalized.moran.mc(crime, nb2listw(COL.nb, style="W"), nsim=99, zero.policy=TRUE,
                    na.action=na.exclude)
#> 
#>  Monte-Carlo simulation of Normalized Moran I
#> 
#> data:  crime 
#> weights: nb2listw(COL.nb, style = "W") 
#> omitted: 5, 7, 9, 15, 19, 33, 40, 45, 47, 48 
#> number of simulations + 1: 100 
#> 
#> statistic = 0.60616, observed rank = 100, p-value = 0.01
#> alternative hypothesis: greater
#> 
normalized.moran.mc(crime, nb2listw(COL.nb, style="W"), nsim=99, zero.policy=TRUE,
                    return_boot=TRUE, na.action=na.exclude)
#> NA observations omitted: 5, 7, 9, 15, 19, 33, 40, 45, 47, 48
#> 
#> DATA PERMUTATION
#> 
#> 
#> Call:
#> boot::boot(data = x, statistic = moran_boot, R = nsim, sim = "permutation", 
#>     listw = listw, zero.policy = zero.policy, parallel = parallel, 
#>     ncpus = ncpus, cl = cl)
#> 
#> 
#> Bootstrap Statistics :
#>      original     bias    std. error
#> t1* 0.6061552 -0.6463009   0.2285857
try(normalized.moran.mc(crime, nb2listw(COL.nb, style="W"), nsim=99, na.action=na.pass))
#> Error in normalized.moran.mc(crime, nb2listw(COL.nb, style = "W"), nsim = 99,  : 
#>   na.pass not permitted
```
