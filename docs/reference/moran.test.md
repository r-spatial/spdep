# Moran's I test for spatial autocorrelation

Moran's test for spatial autocorrelation using a spatial weights matrix
in weights list form. The assumptions underlying the test are sensitive
to the form of the graph of neighbour relationships and other factors,
and results may be checked against those of `moran.mc` permutations.

## Usage

``` r
moran.test(x, listw, randomisation=TRUE, zero.policy=attr(listw, "zero.policy"),
 alternative="greater", rank = FALSE, na.action=na.fail, spChk=NULL,
 adjust.n=TRUE, drop.EI2=FALSE)
```

## Arguments

- x:

  a numeric vector the same length as the neighbours list in listw

- listw:

  a `listw` object created for example by `nb2listw`

- randomisation:

  variance of I calculated under the assumption of randomisation, if
  FALSE normality

- zero.policy:

  default `attr(listw, "zero.policy")` as set when `listw` was created,
  if attribute not set, use global option value; if TRUE assign zero to
  the lagged value of zones without neighbours, if FALSE assign NA

- alternative:

  a character string specifying the alternative hypothesis, must be one
  of greater (default), less or two.sided.

- rank:

  logical value - default FALSE for continuous variables, if TRUE, uses
  the adaptation of Moran's I for ranks suggested by Cliff and Ord
  (1981, p. 46)

- na.action:

  a function (default `na.fail`), can also be `na.omit` or
  `na.exclude` - in these cases the weights list will be subsetted to
  remove NAs in the data. It may be necessary to set zero.policy to TRUE
  because this subsetting may create no-neighbour observations. Note
  that only weights lists created without using the glist argument to
  `nb2listw` may be subsetted. If `na.pass` is used, zero is substituted
  for NA values in calculating the spatial lag

- spChk:

  should the data vector names be checked against the spatial objects
  for identity integrity, TRUE, or FALSE, default NULL to use
  [`get.spChkOption()`](https://r-spatial.github.io/spdep/reference/set.spChkOption.md)

- adjust.n:

  default TRUE, if FALSE the number of observations is not adjusted for
  no-neighbour observations, if TRUE, the number of observations is
  adjusted

- drop.EI2:

  default FALSE, if TRUE, emulate CrimeStat \<= 4.02

## Value

A list with class `htest` containing the following components:

- statistic:

  the value of the standard deviate of Moran's I.

- p.value:

  the p-value of the test.

- estimate:

  the value of the observed Moran's I, its expectation and variance
  under the method assumption.

- alternative:

  a character string describing the alternative hypothesis.

- method:

  a character string giving the assumption used for calculating the
  standard deviate.

- data.name:

  a character string giving the name(s) of the data.

## Note

Var(I) is taken from Cliff and Ord (1969, p. 28), and Goodchild's CATMOG
47 (1986), see also Upton & Fingleton (1985) p. 171; it agrees with
SpaceStat, see Tutorial workbook Chapter 22; VI is the second crude
moment minus the square of the first crude moment. The derivation of the
test (Cliff and Ord, 1981, p. 18) assumes that the weights matrix is
symmetric. For inherently non-symmetric matrices, such as k-nearest
neighbour matrices,
[`listw2U()`](https://r-spatial.github.io/spdep/reference/nb2listw.md)
can be used to make the matrix symmetric.

## References

Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 21; Bivand RS,
Wong DWS 2018 Comparing implementations of global and local indicators
of spatial association. TEST, 27(3), 716–748
[doi:10.1007/s11749-018-0599-x](https://doi.org/10.1007/s11749-018-0599-x)

## Author

Roger Bivand <Roger.Bivand@nhh.no>

## See also

[`moran`](https://r-spatial.github.io/spdep/reference/moran.md),
[`moran.mc`](https://r-spatial.github.io/spdep/reference/moran.mc.md),
[`listw2U`](https://r-spatial.github.io/spdep/reference/nb2listw.md)

## Examples

``` r
data(oldcol)
coords.OLD <- cbind(COL.OLD$X, COL.OLD$Y)
moran.test(COL.OLD$CRIME, nb2listw(COL.nb, style="W"))
#> 
#>  Moran I test under randomisation
#> 
#> data:  COL.OLD$CRIME  
#> weights: nb2listw(COL.nb, style = "W")    
#> 
#> Moran I statistic standard deviate = 5.6341, p-value = 8.797e-09
#> alternative hypothesis: greater
#> sample estimates:
#> Moran I statistic       Expectation          Variance 
#>       0.510951264      -0.020833333       0.008908762 
#> 
moran.test(COL.OLD$CRIME, nb2listw(COL.nb, style="B"))
#> 
#>  Moran I test under randomisation
#> 
#> data:  COL.OLD$CRIME  
#> weights: nb2listw(COL.nb, style = "B")    
#> 
#> Moran I statistic standard deviate = 6.2116, p-value = 2.622e-10
#> alternative hypothesis: greater
#> sample estimates:
#> Moran I statistic       Expectation          Variance 
#>        0.52063815       -0.02083333        0.00759872 
#> 
moran.test(COL.OLD$CRIME, nb2listw(COL.nb, style="C"))
#> 
#>  Moran I test under randomisation
#> 
#> data:  COL.OLD$CRIME  
#> weights: nb2listw(COL.nb, style = "C")    
#> 
#> Moran I statistic standard deviate = 6.2116, p-value = 2.622e-10
#> alternative hypothesis: greater
#> sample estimates:
#> Moran I statistic       Expectation          Variance 
#>        0.52063815       -0.02083333        0.00759872 
#> 
moran.test(COL.OLD$CRIME, nb2listw(COL.nb, style="S"))
#> 
#>  Moran I test under randomisation
#> 
#> data:  COL.OLD$CRIME  
#> weights: nb2listw(COL.nb, style = "S")    
#> 
#> Moran I statistic standard deviate = 5.9786, p-value = 1.125e-09
#> alternative hypothesis: greater
#> sample estimates:
#> Moran I statistic       Expectation          Variance 
#>       0.512957470      -0.020833333       0.007971504 
#> 
moran.test(COL.OLD$CRIME, nb2listw(COL.nb, style="W"),
 randomisation=FALSE)
#> 
#>  Moran I test under normality
#> 
#> data:  COL.OLD$CRIME  
#> weights: nb2listw(COL.nb, style = "W")    
#> 
#> Moran I statistic standard deviate = 5.6754, p-value = 6.92e-09
#> alternative hypothesis: greater
#> sample estimates:
#> Moran I statistic       Expectation          Variance 
#>       0.510951264      -0.020833333       0.008779831 
#> 
colold.lags <- nblag(COL.nb, 3)
moran.test(COL.OLD$CRIME, nb2listw(colold.lags[[2]],
 style="W"))
#> 
#>  Moran I test under randomisation
#> 
#> data:  COL.OLD$CRIME  
#> weights: nb2listw(colold.lags[[2]], style = "W")    
#> 
#> Moran I statistic standard deviate = 2.6076, p-value = 0.004559
#> alternative hypothesis: greater
#> sample estimates:
#> Moran I statistic       Expectation          Variance 
#>       0.168485742      -0.020833333       0.005271314 
#> 
moran.test(COL.OLD$CRIME, nb2listw(colold.lags[[3]],
 style="W"))
#> 
#>  Moran I test under randomisation
#> 
#> data:  COL.OLD$CRIME  
#> weights: nb2listw(colold.lags[[3]], style = "W")    
#> 
#> Moran I statistic standard deviate = -1.7896, p-value = 0.9632
#> alternative hypothesis: greater
#> sample estimates:
#> Moran I statistic       Expectation          Variance 
#>      -0.138930745      -0.020833333       0.004354683 
#> 
print(is.symmetric.nb(COL.nb))
#> [1] TRUE
COL.k4.nb <- knn2nb(knearneigh(coords.OLD, 4))
print(is.symmetric.nb(COL.k4.nb))
#> [1] FALSE
moran.test(COL.OLD$CRIME, nb2listw(COL.k4.nb, style="W"))
#> 
#>  Moran I test under randomisation
#> 
#> data:  COL.OLD$CRIME  
#> weights: nb2listw(COL.k4.nb, style = "W")    
#> 
#> Moran I statistic standard deviate = 7.2183, p-value = 2.632e-13
#> alternative hypothesis: greater
#> sample estimates:
#> Moran I statistic       Expectation          Variance 
#>       0.624933667      -0.020833333       0.008003503 
#> 
moran.test(COL.OLD$CRIME, nb2listw(COL.k4.nb, style="W"),
 randomisation=FALSE)
#> 
#>  Moran I test under normality
#> 
#> data:  COL.OLD$CRIME  
#> weights: nb2listw(COL.k4.nb, style = "W")    
#> 
#> Moran I statistic standard deviate = 7.2711, p-value = 1.782e-13
#> alternative hypothesis: greater
#> sample estimates:
#> Moran I statistic       Expectation          Variance 
#>       0.624933667      -0.020833333       0.007887613 
#> 
cat("Note: non-symmetric weights matrix, use listw2U()")
#> Note: non-symmetric weights matrix, use listw2U()
moran.test(COL.OLD$CRIME, listw2U(nb2listw(COL.k4.nb,
 style="W")))
#> 
#>  Moran I test under randomisation
#> 
#> data:  COL.OLD$CRIME  
#> weights: listw2U(nb2listw(COL.k4.nb, style = "W"))    
#> 
#> Moran I statistic standard deviate = 7.2183, p-value = 2.632e-13
#> alternative hypothesis: greater
#> sample estimates:
#> Moran I statistic       Expectation          Variance 
#>       0.624933667      -0.020833333       0.008003503 
#> 
moran.test(COL.OLD$CRIME, listw2U(nb2listw(COL.k4.nb,
 style="W")), randomisation=FALSE)
#> 
#>  Moran I test under normality
#> 
#> data:  COL.OLD$CRIME  
#> weights: listw2U(nb2listw(COL.k4.nb, style = "W"))    
#> 
#> Moran I statistic standard deviate = 7.2711, p-value = 1.782e-13
#> alternative hypothesis: greater
#> sample estimates:
#> Moran I statistic       Expectation          Variance 
#>       0.624933667      -0.020833333       0.007887613 
#> 
ranks <- rank(COL.OLD$CRIME)
names(ranks) <- rownames(COL.OLD)
moran.test(ranks, nb2listw(COL.nb, style="W"), rank=TRUE)
#> 
#>  Moran I test under randomisation
#> 
#> data:  ranks using rank correction 
#> weights: nb2listw(COL.nb, style = "W")    
#> 
#> Moran I statistic standard deviate = 6.3815, p-value = 8.766e-11
#> alternative hypothesis: greater
#> sample estimates:
#> Moran I statistic       Expectation          Variance 
#>       0.584333495      -0.020833333       0.008992923 
#> 
crime <- COL.OLD$CRIME
is.na(crime) <- sample(1:length(crime), 10)
res <- try(moran.test(crime, nb2listw(COL.nb, style="W"),
 na.action=na.fail))
#> Error in na.fail.default(x) : missing values in object
moran.test(crime, nb2listw(COL.nb, style="W"), zero.policy=TRUE,
 na.action=na.omit)
#> Warning: subsetting caused increase in subgraph count
#> 
#>  Moran I test under randomisation
#> 
#> data:  crime  
#> weights: nb2listw(COL.nb, style = "W") 
#> omitted: 7, 10, 12, 14, 18, 26, 28, 41, 42, 49 
#> n reduced by no-neighbour observations  
#> 
#> Moran I statistic standard deviate = 4.327, p-value = 7.558e-06
#> alternative hypothesis: greater
#> sample estimates:
#> Moran I statistic       Expectation          Variance 
#>        0.47702353       -0.02702703        0.01356988 
#> 
moran.test(crime, nb2listw(COL.nb, style="W"), zero.policy=TRUE,
 na.action=na.exclude)
#> Warning: subsetting caused increase in subgraph count
#> 
#>  Moran I test under randomisation
#> 
#> data:  crime  
#> weights: nb2listw(COL.nb, style = "W") 
#> omitted: 7, 10, 12, 14, 18, 26, 28, 41, 42, 49 
#> n reduced by no-neighbour observations  
#> 
#> Moran I statistic standard deviate = 4.327, p-value = 7.558e-06
#> alternative hypothesis: greater
#> sample estimates:
#> Moran I statistic       Expectation          Variance 
#>        0.47702353       -0.02702703        0.01356988 
#> 
moran.test(crime, nb2listw(COL.nb, style="W"), na.action=na.pass)
#> Warning: NAs in lagged values
#> 
#>  Moran I test under randomisation
#> 
#> data:  crime  
#> weights: nb2listw(COL.nb, style = "W")    
#> 
#> Moran I statistic standard deviate = 2.5935, p-value = 0.00475
#> alternative hypothesis: greater
#> sample estimates:
#> Moran I statistic       Expectation          Variance 
#>       0.222432832      -0.020833333       0.008798241 
#> 
columbus <- st_read(system.file("shapes/columbus.gpkg", package="spData")[1], quiet=TRUE)
col_geoms <- st_geometry(columbus)
col_geoms[1] <- st_buffer(col_geoms[1], dist=-0.05)
st_geometry(columbus) <- col_geoms
(nb1 <- poly2nb(columbus))
#> Warning: some observations have no neighbours;
#> if this seems unexpected, try increasing the snap argument.
#> Warning: neighbour object has 2 sub-graphs;
#> if this sub-graph count seems unexpected, try increasing the snap argument.
#> Neighbour list object:
#> Number of regions: 49 
#> Number of nonzero links: 232 
#> Percentage nonzero weights: 9.662641 
#> Average number of links: 4.734694 
#> 1 region with no links:
#> 1
#> 2 disjoint connected subgraphs
try(lw <- nb2listw(nb1, style="W"))
#> Error in nb2listw(nb1, style = "W") : 
#>   Empty neighbour sets found (zero.policy: FALSE)
(lw <- nb2listw(nb1, style="W", zero.policy=TRUE))
#> Characteristics of weights list object:
#> Neighbour list object:
#> Number of regions: 49 
#> Number of nonzero links: 232 
#> Percentage nonzero weights: 9.662641 
#> Average number of links: 4.734694 
#> 1 region with no links:
#> 1
#> 2 disjoint connected subgraphs
#> 
#> Weights style: W 
#> Weights constants summary:
#>    n   nn S0       S1       S2
#> W 48 2304 48 22.23035 199.8188
moran.test(COL.OLD$CRIME, lw)
#> 
#>  Moran I test under randomisation
#> 
#> data:  COL.OLD$CRIME  
#> weights: lw  
#> n reduced by no-neighbour observations  
#> 
#> Moran I statistic standard deviate = 2.7707, p-value = 0.002797
#> alternative hypothesis: greater
#> sample estimates:
#> Moran I statistic       Expectation          Variance 
#>        0.23901261       -0.02127660        0.00882535 
#> 
```
