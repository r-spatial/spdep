# Permutation test for empirical Bayes index

An empirical Bayes index modification of Moran's I for testing for
spatial autocorrelation in a rate, typically the number of observed
cases in a population at risk. The index value is tested by using nsim
random permutations of the index for the given spatial weighting scheme,
to establish the rank of the observed statistic in relation to the nsim
simulated values.

## Usage

``` r
EBImoran.mc(n, x, listw, nsim, zero.policy = attr(listw, "zero.policy"), 
 alternative = "greater", spChk=NULL, return_boot=FALSE,
 subtract_mean_in_numerator=TRUE)
```

## Arguments

- n:

  a numeric vector of counts of cases the same length as the neighbours
  list in listw

- x:

  a numeric vector of populations at risk the same length as the
  neighbours list in listw

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
  of "greater" (default), "two.sided", or "less"

- spChk:

  should the data vector names be checked against the spatial objects
  for identity integrity, TRUE, or FALSE, default NULL to use
  [`get.spChkOption()`](https://r-spatial.github.io/spdep/reference/set.spChkOption.md)

- return_boot:

  return an object of class `boot` from the equivalent permutation
  bootstrap rather than an object of class `htest`

- subtract_mean_in_numerator:

  default TRUE, if TRUE subtract mean of z in numerator of EBI equation
  on p. 2157 in reference (consulted with Renato Assunção 2016-02-19);
  until February 2016 the default was FALSE agreeing with the printed
  paper.

## Details

The statistic used is (m is the number of observations): \$\$EBI =
\frac{m}{\sum\_{i=1}^{m}\sum\_{j=1}^{m}w\_{ij}}
\frac{\sum\_{i=1}^{m}\sum\_{j=1}^{m}w\_{ij}z_i
z_j}{\sum\_{i=1}^{m}(z_i - \bar{z})^2} \$\$ where: \$\$z_i = \frac{p_i -
b}{\sqrt{v_i}}\$\$ and: \$\$p_i = n_i / x_i\$\$ \$\$v_i = a + (b /
x_i)\$\$ \$\$b = \sum\_{i=1}^{m} n_i / \sum\_{i=1}^{m} x_i \$\$ \$\$a =
s^2 - b / (\sum\_{i=1}^{m} x_i / m)\$\$ \$\$s^2 = \sum\_{i=1}^{m} x_i
(p_i - b)^2 / \sum\_{i=1}^{m} x_i \$\$

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

- z:

  a numerical vector of Empirical Bayes indices as z above

## References

Assunção RM, Reis EA 1999 A new proposal to adjust Moran's I for
population density. Statistics in Medicine 18, pp. 2147–2162; Bivand RS,
Wong DWS 2018 Comparing implementations of global and local indicators
of spatial association. TEST, 27(3), 716–748
[doi:10.1007/s11749-018-0599-x](https://doi.org/10.1007/s11749-018-0599-x)

## Author

Roger Bivand <Roger.Bivand@nhh.no>

## See also

[`moran`](https://r-spatial.github.io/spdep/reference/moran.md),
[`moran.mc`](https://r-spatial.github.io/spdep/reference/moran.mc.md),
[`EBest`](https://r-spatial.github.io/spdep/reference/EBest.md)

## Examples

``` r
nc.sids <- st_read(system.file("shapes/sids.gpkg", package="spData")[1], quiet=TRUE)
rn <- as.character(nc.sids$FIPS)
ncCC89_nb <- read.gal(system.file("weights/ncCC89.gal", package="spData")[1],
 region.id=rn)
#> Warning: neighbour object has 3 sub-graphs
EBImoran.mc(nc.sids$SID74, nc.sids$BIR74,
 nb2listw(ncCC89_nb, style="B", zero.policy=TRUE), nsim=999,
 alternative="two.sided", zero.policy=TRUE)
#> The default for subtract_mean_in_numerator set TRUE from February 2016
#> 
#>  Monte-Carlo simulation of Empirical Bayes Index (mean subtracted)
#> 
#> data:  cases: nc.sids$SID74, risk population: nc.sids$BIR74
#> weights: nb2listw(ncCC89_nb, style = "B", zero.policy = TRUE)
#> number of simulations + 1: 1000
#> 
#> statistic = 0.25789, observed rank = 999, p-value = 0.002
#> alternative hypothesis: two.sided
#> 
sids.p <- nc.sids$SID74 / nc.sids$BIR74
moran.mc(sids.p, nb2listw(ncCC89_nb, style="B", zero.policy=TRUE),
 nsim=999, alternative="two.sided", zero.policy=TRUE)
#> 
#>  Monte-Carlo simulation of Moran I
#> 
#> data:  sids.p 
#> weights: nb2listw(ncCC89_nb, style = "B", zero.policy = TRUE)  
#> number of simulations + 1: 1000 
#> 
#> statistic = 0.20904, observed rank = 997, p-value = 0.006
#> alternative hypothesis: two.sided
#> 
```
