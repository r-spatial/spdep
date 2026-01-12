# Permutation test for Lee's L statistic

A permutation test for Lee's L statistic calculated by using nsim random
permutations of x and y for the given spatial weighting scheme, to
establish the rank of the observed statistic in relation to the nsim
simulated values.

## Usage

``` r
lee.mc(x, y, listw, nsim, zero.policy=attr(listw, "zero.policy"), alternative="greater",
 na.action=na.fail, spChk=NULL, return_boot=FALSE)
```

## Arguments

- x:

  a numeric vector the same length as the neighbours list in listw

- y:

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

## Value

A list with class `htest` and `mc.sim` containing the following
components:

- statistic:

  the value of the observed Lee's L.

- parameter:

  the rank of the observed Lee's L.

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

Lee (2001). Developing a bivariate spatial association measure: An
integration of Pearson's r and Moran's I. J Geograph Syst 3: 369-385

## Author

Roger Bivand, Virgilio GÃ³mez-Rubio <Virgilio.Gomez@uclm.es>

## See also

[`lee`](https://r-spatial.github.io/spdep/reference/lee.md)

## Examples

``` r
data(boston, package="spData")
lw<-nb2listw(boston.soi)

x<-boston.c$CMEDV
y<-boston.c$CRIM

lee.mc(x, y, nsim=99, lw, zero.policy=TRUE, alternative="two.sided")
#> 
#>  Monte-Carlo simulation of Lee's L
#> 
#> data:  x ,  y 
#> weights: lw  
#> number of simulations + 1: 100 
#> 
#> statistic = -0.3263, observed rank = 1, p-value = 0.02
#> alternative hypothesis: two.sided
#> 

#Test with missing values
x[1:5]<-NA
y[3:7]<-NA

lee.mc(x, y, nsim=99, lw, zero.policy=TRUE, alternative="two.sided", 
   na.action=na.omit)
#> 
#>  Monte-Carlo simulation of Lee's L
#> 
#> data:  x ,  y 
#> weights: lw 
#> omitted: 1, 2, 3, 4, 5, 6, 7 
#> number of simulations + 1: 100 
#> 
#> statistic = -0.32447, observed rank = 1, p-value = 0.02
#> alternative hypothesis: two.sided
#> 
```
