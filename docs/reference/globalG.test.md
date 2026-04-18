# Global G test for spatial autocorrelation

The global G statistic for spatial autocorrelation, complementing the
local Gi LISA measures:
[`localG`](https://r-spatial.github.io/spdep/reference/localG.md).

## Usage

``` r
globalG.test(x, listw, zero.policy=attr(listw, "zero.policy"), alternative="greater",
 spChk=NULL, adjust.n=TRUE, B1correct=TRUE, adjust.x=TRUE, Arc_all_x=FALSE,
 na.action=na.fail)
```

## Arguments

- x:

  a numeric vector the same length as the neighbours list in listw

- listw:

  a `listw` object created for example by `nb2listw`; if a sequence of
  distance bands is to be used, it is recommended that the weights style
  be binary (one of `c("B", "C", "U")`).

- zero.policy:

  default `attr(listw, "zero.policy")` as set when `listw` was created,
  if attribute not set, use global option value; if TRUE assign zero to
  the lagged value of zones without neighbours, if FALSE assign NA

- alternative:

  a character string specifying the alternative hypothesis, must be one
  of "greater" (default), "less" or "two.sided".

- spChk:

  should the data vector names be checked against the spatial objects
  for identity integrity, TRUE, or FALSE, default NULL to use
  [`get.spChkOption()`](https://r-spatial.github.io/spdep/reference/set.spChkOption.md)

- adjust.n:

  default TRUE, if FALSE the number of observations is not adjusted for
  no-neighbour observations, if TRUE, the number of observations is
  adjusted

- B1correct:

  default TRUE, if TRUE, the erratum referenced below: "On page 195, the
  coefficient of W2 in B1, (just below center of the page) should be 6,
  not 3." is applied; if FALSE, 3 is used (as in CrimeStat IV)

- adjust.x:

  default TRUE, if TRUE, x values of observations with no neighbours are
  omitted in the denominator of G

- Arc_all_x:

  default FALSE, if Arc_all_x=TRUE and adjust.x=TRUE, use the full x
  vector in part of the denominator term for G

- na.action:

  a function (default `na.fail`), can also be `na.omit` or
  `na.exclude` - in these cases the weights list will be subsetted to
  remove NAs in the data. It may be necessary to set zero.policy to TRUE
  because this subsetting may create no-neighbour observations. Note
  that only weights lists created without using the glist argument to
  `nb2listw` may be subsetted. `na.pass` is not permitted.

## Value

A list with class `htest` containing the following components:

- statistic:

  the value of the standard deviate of Getis-Ord G.

- p.value:

  the p-value of the test.

- estimate:

  the value of the observed statistic, its expectation and variance.

- alternative:

  a character string describing the alternative hypothesis.

- data.name:

  a character string giving the name(s) of the data.

## References

Getis. A, Ord, J. K. 1992 The analysis of spatial association by use of
distance statistics, *Geographical Analysis*, 24, p. 195; see also
Getis. A, Ord, J. K. 1993 Erratum, *Geographical Analysis*, 25, p. 276;
Bivand RS, Wong DWS 2018 Comparing implementations of global and local
indicators of spatial association. TEST, 27(3), 716–748
[doi:10.1007/s11749-018-0599-x](https://doi.org/10.1007/s11749-018-0599-x)

## Author

Hisaji ONO <hi-ono@mn.xdsl.ne.jp> and Roger Bivand <Roger.Bivand@nhh.no>

## See also

[`localG`](https://r-spatial.github.io/spdep/reference/localG.md)

## Examples

``` r
nc.sids <- st_read(system.file("shapes/sids.gpkg", package="spData")[1], quiet=TRUE)
sidsrate79 <- (1000*nc.sids$SID79)/nc.sids$BIR79
dists <- c(10, 20, 30, 33, 40, 50, 60, 70, 80, 90, 100)
ndists <- length(dists)
ZG <- vector(mode="list", length=ndists)
names(ZG) <- as.character(dists)
milesxy <- cbind(nc.sids$east, nc.sids$north)
for (i in 1:ndists) {
  thisnb <- dnearneigh(milesxy, 0, dists[i])
  thislw <- nb2listw(thisnb, style="B", zero.policy=TRUE)
  ZG[[i]] <- globalG.test(sidsrate79, thislw, zero.policy=TRUE)
}
#> Warning: neighbour object has 98 sub-graphs
#> Warning: neighbour object has 37 sub-graphs
#> Warning: neighbour object has 3 sub-graphs
t(sapply(ZG, function(x) c(x$estimate[1], x$statistic, p.value=unname(x$p.value))))
#>     Global G statistic standard deviate   p.value
#> 10          0.33548581       0.04859237 0.4806221
#> 20          0.02024945      -0.73262399 0.7681061
#> 30          0.04032434      -0.75011405 0.7734070
#> 33          0.05312271       0.40157023 0.3440002
#> 40          0.07400279      -0.04345713 0.5173314
#> 50          0.11471743       0.58686472 0.2786473
#> 60          0.15457553      -0.35823892 0.6399177
#> 70          0.19839023      -0.27864299 0.6097406
#> 80          0.24606972      -0.18791364 0.5745278
#> 90          0.30073463       0.11457610 0.4543906
#> 100         0.34879996       0.31591356 0.3760341
for (i in 1:ndists) {
  thisnb <- dnearneigh(milesxy, 0, dists[i])
  thislw <- nb2listw(thisnb, style="B", zero.policy=TRUE)
  ZG[[i]] <- globalG.test(sidsrate79, thislw, zero.policy=TRUE, alternative="two.sided")
}
#> Warning: neighbour object has 98 sub-graphs
#> Warning: neighbour object has 37 sub-graphs
#> Warning: neighbour object has 3 sub-graphs
t(sapply(ZG, function(x) c(x$estimate[1], x$statistic, p.value=unname(x$p.value))))
#>     Global G statistic standard deviate   p.value
#> 10          0.33548581       0.04859237 0.9612441
#> 20          0.02024945      -0.73262399 0.4637878
#> 30          0.04032434      -0.75011405 0.4531860
#> 33          0.05312271       0.40157023 0.6880003
#> 40          0.07400279      -0.04345713 0.9653371
#> 50          0.11471743       0.58686472 0.5572946
#> 60          0.15457553      -0.35823892 0.7201645
#> 70          0.19839023      -0.27864299 0.7805188
#> 80          0.24606972      -0.18791364 0.8509443
#> 90          0.30073463       0.11457610 0.9087811
#> 100         0.34879996       0.31591356 0.7520681
data(oldcol)
crime <- COL.OLD$CRIME
is.na(crime) <- sample(1:length(crime), 10)
res <- try(globalG.test(crime, nb2listw(COL.nb, style="B"),
 na.action=na.fail))
#> Error in na.fail.default(x) : missing values in object
globalG.test(crime, nb2listw(COL.nb, style="B"), zero.policy=TRUE,
 na.action=na.omit)
#> Warning: subsetting caused increase in subgraph count
#> 
#>  Getis-Ord global G statistic
#> 
#> data:  crime 
#> weights: nb2listw(COL.nb, style = "B") 
#> omitted: 2, 3, 5, 6, 7, 10, 12, 30, 32, 42 
#> n reduced by no-neighbour observations 
#> 
#> standard deviate = 2.2378, p-value = 0.01262
#> alternative hypothesis: greater
#> sample estimates:
#> Global G statistic        Expectation           Variance 
#>       1.207385e-01       1.038407e-01       5.701826e-05 
#> 
globalG.test(crime, nb2listw(COL.nb, style="B"), zero.policy=TRUE,
 na.action=na.exclude)
#> Warning: subsetting caused increase in subgraph count
#> 
#>  Getis-Ord global G statistic
#> 
#> data:  crime 
#> weights: nb2listw(COL.nb, style = "B") 
#> omitted: 2, 3, 5, 6, 7, 10, 12, 30, 32, 42 
#> n reduced by no-neighbour observations 
#> 
#> standard deviate = 2.2378, p-value = 0.01262
#> alternative hypothesis: greater
#> sample estimates:
#> Global G statistic        Expectation           Variance 
#>       1.207385e-01       1.038407e-01       5.701826e-05 
#> 
try(globalG.test(crime, nb2listw(COL.nb, style="B"), na.action=na.pass))
#> Error in globalG.test(crime, nb2listw(COL.nb, style = "B"), na.action = na.pass) : 
#>   na.pass not permitted
```
