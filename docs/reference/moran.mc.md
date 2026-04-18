# Permutation test for Moran's I statistic

A permutation test for Moran's I statistic calculated by using nsim
random permutations of x for the given spatial weighting scheme, to
establish the rank of the observed statistic in relation to the nsim
simulated values.

## Usage

``` r
moran.mc(x, listw, nsim, zero.policy=attr(listw, "zero.policy"),
 alternative="greater", na.action=na.fail, spChk=NULL, return_boot=FALSE,
 adjust.n=TRUE)
```

## Arguments

- x:

  a numeric vector the same length as the neighbours list in listw

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

- na.action:

  a function (default `na.fail`), can also be `na.omit` or
  `na.exclude` - in these cases the weights list will be subsetted to
  remove NAs in the data. It may be necessary to set zero.policy to TRUE
  because this subsetting may create no-neighbour observations. Note
  that only weights lists created without using the glist argument to
  `nb2listw` may be subsetted. `na.pass` is not permitted because it is
  meaningless in a permutation test.

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

- statistic:

  the value of the observed Moran's I.

- parameter:

  the rank of the observed Moran's I.

- p.value:

  the pseudo p-value of the test.

- alternative:

  a character string describing the alternative hypothesis.

- method:

  a character string giving the method used.

- data.name:

  a character string giving the name(s) of the data, and the number of
  simulations.

- res:

  nsim simulated values of statistic, final value is observed statistic

## References

Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 63-5.

## Author

Roger Bivand <Roger.Bivand@nhh.no>

## See also

[`moran`](https://r-spatial.github.io/spdep/reference/moran.md),
[`moran.test`](https://r-spatial.github.io/spdep/reference/moran.test.md)

## Examples

``` r
data(oldcol)
colw <- nb2listw(COL.nb, style="W")
nsim <- 99
set.seed(1234)
sim1 <- moran.mc(COL.OLD$CRIME, listw=colw, nsim=nsim)
sim1
#> 
#>  Monte-Carlo simulation of Moran I
#> 
#> data:  COL.OLD$CRIME 
#> weights: colw  
#> number of simulations + 1: 100 
#> 
#> statistic = 0.51095, observed rank = 100, p-value = 0.01
#> alternative hypothesis: greater
#> 
mean(sim1$res[1:nsim])
#> [1] -0.003196117
var(sim1$res[1:nsim])
#> [1] 0.008543554
summary(sim1$res[1:nsim])
#>      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#> -0.180085 -0.073721 -0.009528 -0.003196  0.063865  0.206779 
colold.lags <- nblag(COL.nb, 3)
set.seed(1234)
sim2 <- moran.mc(COL.OLD$CRIME, nb2listw(colold.lags[[2]],
 style="W"), nsim=nsim)
summary(sim2$res[1:nsim])
#>     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#> -0.15905 -0.07503 -0.01535 -0.02060  0.03305  0.13023 
sim3 <- moran.mc(COL.OLD$CRIME, nb2listw(colold.lags[[3]],
 style="W"), nsim=nsim)
summary(sim3$res[1:nsim])
#>      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#> -0.192986 -0.055468 -0.024154 -0.026327  0.004769  0.129996 
crime <- COL.OLD$CRIME
is.na(crime) <- sample(1:length(crime), 10)
try(moran.mc(crime, nb2listw(COL.nb, style="W"), nsim=99,
 na.action=na.fail))
#> Error in na.fail.default(x) : missing values in object
moran.mc(crime, nb2listw(COL.nb, style="W"), nsim=99, zero.policy=TRUE,
 na.action=na.omit)
#> 
#>  Monte-Carlo simulation of Moran I
#> 
#> data:  crime 
#> weights: nb2listw(COL.nb, style = "W") 
#> omitted: 5, 7, 9, 15, 19, 33, 40, 45, 47, 48 
#> number of simulations + 1: 100 
#> 
#> statistic = 0.42342, observed rank = 100, p-value = 0.01
#> alternative hypothesis: greater
#> 
moran.mc(crime, nb2listw(COL.nb, style="W"), nsim=99, zero.policy=TRUE,
 return_boot=TRUE, na.action=na.omit)
#> NA observations omitted: 5, 7, 9, 15, 19, 33, 40, 45, 47, 48
#> 
#> DATA PERMUTATION
#> 
#> 
#> Call:
#> boot(data = x, statistic = moran_boot, R = nsim, sim = "permutation", 
#>     listw = listw, n = n, S0 = S0, zero.policy = zero.policy, 
#>     parallel = parallel, ncpus = ncpus, cl = cl)
#> 
#> 
#> Bootstrap Statistics :
#>      original     bias    std. error
#> t1* 0.4234221 -0.4429891   0.1293935
moran.mc(crime, nb2listw(COL.nb, style="W"), nsim=99, zero.policy=TRUE,
 na.action=na.exclude)
#> 
#>  Monte-Carlo simulation of Moran I
#> 
#> data:  crime 
#> weights: nb2listw(COL.nb, style = "W") 
#> omitted: 5, 7, 9, 15, 19, 33, 40, 45, 47, 48 
#> number of simulations + 1: 100 
#> 
#> statistic = 0.42342, observed rank = 100, p-value = 0.01
#> alternative hypothesis: greater
#> 
moran.mc(crime, nb2listw(COL.nb, style="W"), nsim=99, zero.policy=TRUE,
 return_boot=TRUE, na.action=na.exclude)
#> NA observations omitted: 5, 7, 9, 15, 19, 33, 40, 45, 47, 48
#> 
#> DATA PERMUTATION
#> 
#> 
#> Call:
#> boot(data = x, statistic = moran_boot, R = nsim, sim = "permutation", 
#>     listw = listw, n = n, S0 = S0, zero.policy = zero.policy, 
#>     parallel = parallel, ncpus = ncpus, cl = cl)
#> 
#> 
#> Bootstrap Statistics :
#>      original     bias    std. error
#> t1* 0.4234221 -0.4384413   0.1243937
try(moran.mc(crime, nb2listw(COL.nb, style="W"), nsim=99, na.action=na.pass))
#> Error in moran.mc(crime, nb2listw(COL.nb, style = "W"), nsim = 99, na.action = na.pass) : 
#>   na.pass not permitted
```
