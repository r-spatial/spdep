# Permutation test for Geary's C statistic

A permutation test for Geary's C statistic calculated by using nsim
random permutations of x for the given spatial weighting scheme, to
establish the rank of the observed statistic in relation to the nsim
simulated values.

## Usage

``` r
geary.mc(x, listw, nsim, zero.policy=attr(listw, "zero.policy"), alternative="greater",
 spChk=NULL, adjust.n=TRUE, return_boot=FALSE, na.action=na.fail, scale=TRUE)
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
  of "greater" (default), or "less"; this reversal corresponds to that
  on
  [`geary.test`](https://r-spatial.github.io/spdep/reference/geary.test.md)
  described in the section on the output statistic value, based on Cliff
  and Ord 1973, p. 21 (changed 2011-04-11, thanks to Daniel Garavito).

- spChk:

  should the data vector names be checked against the spatial objects
  for identity integrity, TRUE, or FALSE, default NULL to use
  [`get.spChkOption()`](https://r-spatial.github.io/spdep/reference/set.spChkOption.md)

- adjust.n:

  default TRUE, if FALSE the number of observations is not adjusted for
  no-neighbour observations, if TRUE, the number of observations is
  adjusted

- return_boot:

  return an object of class `boot` from the equivalent permutation
  bootstrap rather than an object of class `htest`

- na.action:

  a function (default `na.fail`), can also be `na.omit` or
  `na.exclude` - in these cases the weights list will be subsetted to
  remove NAs in the data. It may be necessary to set zero.policy to TRUE
  because this subsetting may create no-neighbour observations. Note
  that only weights lists created without using the glist argument to
  `nb2listw` may be subsetted. `na.pass` is not permitted because it is
  meaningless in a permutation test.

- scale:

  default TRUE, may be FALSE to revert changes made to accommodate
  `localC` in November 2021 (see \#151)

## Value

A list with class `htest` and `mc.sim` containing the following
components:

- statistic:

  the value of the observed Geary's C.

- parameter:

  the rank of the observed Geary's C.

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

[`geary`](https://r-spatial.github.io/spdep/reference/geary.md),
[`geary.test`](https://r-spatial.github.io/spdep/reference/geary.test.md)

## Examples

``` r
data(oldcol)
set.seed(1)
sim1 <- geary.mc(COL.OLD$CRIME, nb2listw(COL.nb, style="W"),
 nsim=99, alternative="less")
sim1
#> 
#>  Monte-Carlo simulation of Geary C
#> 
#> data:  COL.OLD$CRIME 
#> weights: nb2listw(COL.nb, style = "W")  
#> number of simulations + 1: 100 
#> 
#> statistic = 0.52987, observed rank = 1, p-value = 0.99
#> alternative hypothesis: less
#> 
mean(sim1$res)
#> [1] 1.002705
var(sim1$res)
#> [1] 0.01043616
summary(sim1$res)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.5299  0.9362  1.0141  1.0027  1.0673  1.2459 
colold.lags <- nblag(COL.nb, 3)
sim2 <- geary.mc(COL.OLD$CRIME, nb2listw(colold.lags[[2]],
 style="W"), nsim=99)
sim2
#> 
#>  Monte-Carlo simulation of Geary C
#> 
#> data:  COL.OLD$CRIME 
#> weights: nb2listw(colold.lags[[2]], style = "W")  
#> number of simulations + 1: 100 
#> 
#> statistic = 0.81129, observed rank = 1, p-value = 0.01
#> alternative hypothesis: greater
#> 
summary(sim2$res)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.8113  0.9500  1.0173  1.0147  1.0731  1.1962 
sim3 <- geary.mc(COL.OLD$CRIME, nb2listw(colold.lags[[3]],
 style="W"), nsim=99)
sim3
#> 
#>  Monte-Carlo simulation of Geary C
#> 
#> data:  COL.OLD$CRIME 
#> weights: nb2listw(colold.lags[[3]], style = "W")  
#> number of simulations + 1: 100 
#> 
#> statistic = 1.1303, observed rank = 89, p-value = 0.89
#> alternative hypothesis: greater
#> 
summary(sim3$res)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.7903  0.9609  0.9989  1.0148  1.0656  1.2726 
crime <- COL.OLD$CRIME
is.na(crime) <- sample(1:length(crime), 10)
try(geary.mc(crime, nb2listw(COL.nb, style="W"), nsim=99,
 na.action=na.fail))
#> Error in na.fail.default(x) : missing values in object
geary.mc(crime, nb2listw(COL.nb, style="W"), nsim=99, zero.policy=TRUE,
 na.action=na.omit)
#> Warning: subsetting caused increase in subgraph count
#> 
#>  Monte-Carlo simulation of Geary C
#> 
#> data:  crime 
#> weights: nb2listw(COL.nb, style = "W") 
#> omitted: 14, 16, 17, 21, 24, 31, 32, 42, 44, 47 
#> number of simulations + 1: 100 
#> 
#> statistic = 0.49368, observed rank = 1, p-value = 0.01
#> alternative hypothesis: greater
#> 
geary.mc(crime, nb2listw(COL.nb, style="W"), nsim=99, zero.policy=TRUE,
 return_boot=TRUE, na.action=na.omit)
#> Warning: subsetting caused increase in subgraph count
#> NA observations omitted: 14, 16, 17, 21, 24, 31, 32, 42, 44, 47
#> 
#> DATA PERMUTATION
#> 
#> 
#> Call:
#> boot(data = x, statistic = geary_boot, R = nsim, sim = "permutation", 
#>     listw = listw, n = n, n1 = wc$n1, S0 = wc$S0, zero.policy = zero.policy, 
#>     scale = scale, parallel = parallel, ncpus = ncpus, cl = cl)
#> 
#> 
#> Bootstrap Statistics :
#>      original    bias    std. error
#> t1* 0.4936767 0.4891213   0.1517218
geary.mc(crime, nb2listw(COL.nb, style="W"), nsim=99, zero.policy=TRUE,
 na.action=na.exclude)
#> Warning: subsetting caused increase in subgraph count
#> 
#>  Monte-Carlo simulation of Geary C
#> 
#> data:  crime 
#> weights: nb2listw(COL.nb, style = "W") 
#> omitted: 14, 16, 17, 21, 24, 31, 32, 42, 44, 47 
#> number of simulations + 1: 100 
#> 
#> statistic = 0.49368, observed rank = 1, p-value = 0.01
#> alternative hypothesis: greater
#> 
geary.mc(crime, nb2listw(COL.nb, style="W"), nsim=99, zero.policy=TRUE,
 return_boot=TRUE, na.action=na.exclude)
#> Warning: subsetting caused increase in subgraph count
#> NA observations omitted: 14, 16, 17, 21, 24, 31, 32, 42, 44, 47
#> 
#> DATA PERMUTATION
#> 
#> 
#> Call:
#> boot(data = x, statistic = geary_boot, R = nsim, sim = "permutation", 
#>     listw = listw, n = n, n1 = wc$n1, S0 = wc$S0, zero.policy = zero.policy, 
#>     scale = scale, parallel = parallel, ncpus = ncpus, cl = cl)
#> 
#> 
#> Bootstrap Statistics :
#>      original    bias    std. error
#> t1* 0.4936767 0.4570379   0.1307276
try(geary.mc(crime, nb2listw(COL.nb, style="W"), nsim=99, na.action=na.pass))
#> Error in geary.mc(crime, nb2listw(COL.nb, style = "W"), nsim = 99, na.action = na.pass) : 
#>   na.pass not permitted
```
