# BB join count statistic for k-coloured factors

The BB join count test for spatial autocorrelation using a spatial
weights matrix in weights list form for testing whether same-colour
joins occur more frequently than would be expected if the zones were
labelled in a spatially random way. The assumptions underlying the test
are sensitive to the form of the graph of neighbour relationships and
other factors, and results may be checked against those of
`joincount.mc` permutations.

## Usage

``` r
joincount.test(fx, listw, zero.policy=attr(listw, "zero.policy"), alternative="greater",
 sampling="nonfree", spChk=NULL, adjust.n=TRUE)
# S3 method for class 'jclist'
print(x, ...)
```

## Arguments

- fx:

  a factor of the same length as the neighbours and weights objects in
  listw; use of an ordered factor is not well understood

- listw:

  a `listw` object created for example by `nb2listw`

- zero.policy:

  default `attr(listw, "zero.policy")` as set when `listw` was created,
  if attribute not set, use global option value; if TRUE assign zero to
  the lagged value of zones without neighbours, if FALSE assign NA

- alternative:

  a character string specifying the alternative hypothesis, must be one
  of "greater" (default), "less" or "two.sided".

- sampling:

  default “nonfree”, may be “free”

- adjust.n:

  default TRUE, if FALSE the number of observations is not adjusted for
  no-neighbour observations, if TRUE, the number of observations is
  adjusted consistently (up to and including spdep 0.3-28 the adjustment
  was inconsistent - thanks to Tomoki NAKAYA for a careful bug report)

- spChk:

  should the data vector names be checked against the spatial objects
  for identity integrity, TRUE, or FALSE, default NULL to use
  [`get.spChkOption()`](https://r-spatial.github.io/spdep/reference/set.spChkOption.md)

- x:

  object to be printed

- ...:

  arguments to be passed through for printing

## Value

A list with class `jclist` of lists with class `htest` for each of the k
colours containing the following components:

- statistic:

  the value of the standard deviate of the join count statistic.

- p.value:

  the p-value of the test.

- estimate:

  the value of the observed statistic, its expectation and variance
  under non-free sampling.

- alternative:

  a character string describing the alternative hypothesis.

- method:

  a character string giving the method used.

- data.name:

  a character string giving the name(s) of the data.

## References

Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, pp. 19-20.

## Author

Roger Bivand <Roger.Bivand@nhh.no>

## Note

The derivation of the test (Cliff and Ord, 1981, p. 18) assumes that the
weights matrix is symmetric. For inherently non-symmetric matrices, such
as k-nearest neighbour matrices,
[`listw2U()`](https://r-spatial.github.io/spdep/reference/nb2listw.md)
can be used to make the matrix symmetric. In non-symmetric weights
matrix cases, the variance of the test statistic may be negative.

## See also

[`joincount.mc`](https://r-spatial.github.io/spdep/reference/joincount.mc.md),
[`joincount.multi`](https://r-spatial.github.io/spdep/reference/joincount.multi.md),
[`listw2U`](https://r-spatial.github.io/spdep/reference/nb2listw.md)

## Examples

``` r
data(oldcol)
HICRIME <- cut(COL.OLD$CRIME, breaks=c(0,35,80), labels=c("low","high"))
names(HICRIME) <- rownames(COL.OLD)
joincount.test(HICRIME, nb2listw(COL.nb, style="B"))
#> 
#>  Join count test under nonfree sampling
#> 
#> data:  HICRIME 
#> weights: nb2listw(COL.nb, style = "B") 
#> 
#> Std. deviate for low = 1.0141, p-value = 0.1553
#> alternative hypothesis: greater
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
#> Std. deviate for high = 6.3307, p-value = 1.22e-10
#> alternative hypothesis: greater
#> sample estimates:
#> Same colour statistic           Expectation              Variance 
#>              54.00000              27.22449              17.88838 
#> 
joincount.test(HICRIME, nb2listw(COL.nb, style="B"), sampling="free")
#> 
#>  Join count test under free sampling
#> 
#> data:  HICRIME 
#> weights: nb2listw(COL.nb, style = "B") 
#> 
#> Std. deviate for low = 0.3993, p-value = 0.3448
#> alternative hypothesis: greater
#> sample estimates:
#> Same colour statistic           Expectation              Variance 
#>              34.00000              30.19575              90.76809 
#> 
#> 
#>  Join count test under free sampling
#> 
#> data:  HICRIME 
#> weights: nb2listw(COL.nb, style = "B") 
#> 
#> Std. deviate for high = 2.8518, p-value = 0.002173
#> alternative hypothesis: greater
#> sample estimates:
#> Same colour statistic           Expectation              Variance 
#>               54.0000               27.8284               84.2198 
#> 
joincount.test(HICRIME, nb2listw(COL.nb, style="C"))
#> 
#>  Join count test under nonfree sampling
#> 
#> data:  HICRIME 
#> weights: nb2listw(COL.nb, style = "C") 
#> 
#> Std. deviate for low = 1.0141, p-value = 0.1553
#> alternative hypothesis: greater
#> sample estimates:
#> Same colour statistic           Expectation              Variance 
#>             7.1810345             6.2500000             0.8428969 
#> 
#> 
#>  Join count test under nonfree sampling
#> 
#> data:  HICRIME 
#> weights: nb2listw(COL.nb, style = "C") 
#> 
#> Std. deviate for high = 6.3307, p-value = 1.22e-10
#> alternative hypothesis: greater
#> sample estimates:
#> Same colour statistic           Expectation              Variance 
#>            11.4051724             5.7500000             0.7979712 
#> 
joincount.test(HICRIME, nb2listw(COL.nb, style="S"))
#> 
#>  Join count test under nonfree sampling
#> 
#> data:  HICRIME 
#> weights: nb2listw(COL.nb, style = "S") 
#> 
#> Std. deviate for low = 2.5786, p-value = 0.00496
#> alternative hypothesis: greater
#> sample estimates:
#> Same colour statistic           Expectation              Variance 
#>             8.2425673             6.2500000             0.5971141 
#> 
#> 
#>  Join count test under nonfree sampling
#> 
#> data:  HICRIME 
#> weights: nb2listw(COL.nb, style = "S") 
#> 
#> Std. deviate for high = 6.1736, p-value = 3.337e-10
#> alternative hypothesis: greater
#> sample estimates:
#> Same colour statistic           Expectation              Variance 
#>            10.4249914             5.7500000             0.5734265 
#> 
joincount.test(HICRIME, nb2listw(COL.nb, style="W"))
#> 
#>  Join count test under nonfree sampling
#> 
#> data:  HICRIME 
#> weights: nb2listw(COL.nb, style = "W") 
#> 
#> Std. deviate for low = 4.6675, p-value = 1.524e-06
#> alternative hypothesis: greater
#> sample estimates:
#> Same colour statistic           Expectation              Variance 
#>             9.5190476             6.2500000             0.4905378 
#> 
#> 
#>  Join count test under nonfree sampling
#> 
#> data:  HICRIME 
#> weights: nb2listw(COL.nb, style = "W") 
#> 
#> Std. deviate for high = 5.1205, p-value = 1.523e-07
#> alternative hypothesis: greater
#> sample estimates:
#> Same colour statistic           Expectation              Variance 
#>             9.2920635             5.7500000             0.4784979 
#> 
by(card(COL.nb), HICRIME, summary)
#> HICRIME: low
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>    2.00    2.00    4.00    3.84    4.00   10.00 
#> ------------------------------------------------------------ 
#> HICRIME: high
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   3.000   4.750   6.000   5.667   7.000   9.000 
print(is.symmetric.nb(COL.nb))
#> [1] TRUE
coords.OLD <- cbind(COL.OLD$X, COL.OLD$Y)
COL.k4.nb <- knn2nb(knearneigh(coords.OLD, 4))
print(is.symmetric.nb(COL.k4.nb))
#> [1] FALSE
joincount.test(HICRIME, nb2listw(COL.k4.nb, style="B"))
#> 
#>  Join count test under nonfree sampling
#> 
#> data:  HICRIME 
#> weights: nb2listw(COL.k4.nb, style = "B") 
#> 
#> Std. deviate for low = 4.3698, p-value = 6.217e-06
#> alternative hypothesis: greater
#> sample estimates:
#> Same colour statistic           Expectation              Variance 
#>             36.500000             25.000000              6.925749 
#> 
#> 
#>  Join count test under nonfree sampling
#> 
#> data:  HICRIME 
#> weights: nb2listw(COL.k4.nb, style = "B") 
#> 
#> Std. deviate for high = 6.7293, p-value = 8.523e-12
#> alternative hypothesis: greater
#> sample estimates:
#> Same colour statistic           Expectation              Variance 
#>             40.500000             23.000000              6.762918 
#> 
cat("Note non-symmetric weights matrix - use listw2U()\n")
#> Note non-symmetric weights matrix - use listw2U()
joincount.test(HICRIME, listw2U(nb2listw(COL.k4.nb, style="B")))
#> 
#>  Join count test under nonfree sampling
#> 
#> data:  HICRIME 
#> weights: listw2U(nb2listw(COL.k4.nb, style = "B")) 
#> 
#> Std. deviate for low = 4.3698, p-value = 6.217e-06
#> alternative hypothesis: greater
#> sample estimates:
#> Same colour statistic           Expectation              Variance 
#>             36.500000             25.000000              6.925749 
#> 
#> 
#>  Join count test under nonfree sampling
#> 
#> data:  HICRIME 
#> weights: listw2U(nb2listw(COL.k4.nb, style = "B")) 
#> 
#> Std. deviate for high = 6.7293, p-value = 8.523e-12
#> alternative hypothesis: greater
#> sample estimates:
#> Same colour statistic           Expectation              Variance 
#>             40.500000             23.000000              6.762918 
#> 
columbus <- st_read(system.file("shapes/columbus.gpkg", package="spData")[1], quiet=TRUE)
HICRIME <- cut(columbus$CRIME, breaks=c(0,35,80), labels=c("low","high"))
(nb <- poly2nb(columbus))
#> Neighbour list object:
#> Number of regions: 49 
#> Number of nonzero links: 236 
#> Percentage nonzero weights: 9.829238 
#> Average number of links: 4.816327 
lw <- nb2listw(nb, style="B")
joincount.test(HICRIME, lw)
#> 
#>  Join count test under nonfree sampling
#> 
#> data:  HICRIME 
#> weights: lw 
#> 
#> Std. deviate for low = 1.1164, p-value = 0.1321
#> alternative hypothesis: greater
#> sample estimates:
#> Same colour statistic           Expectation              Variance 
#>              35.00000              30.10204              19.24708 
#> 
#> 
#>  Join count test under nonfree sampling
#> 
#> data:  HICRIME 
#> weights: lw 
#> 
#> Std. deviate for high = 6.163, p-value = 3.57e-10
#> alternative hypothesis: greater
#> sample estimates:
#> Same colour statistic           Expectation              Variance 
#>              54.00000              27.69388              18.21936 
#> 
col_geoms <- st_geometry(columbus)
col_geoms[21] <- st_buffer(col_geoms[21], dist=-0.05)
st_geometry(columbus) <- col_geoms
(nb <- poly2nb(columbus))
#> Warning: some observations have no neighbours;
#> if this seems unexpected, try increasing the snap argument.
#> Warning: neighbour object has 3 sub-graphs;
#> if this sub-graph count seems unexpected, try increasing the snap argument.
#> Neighbour list object:
#> Number of regions: 49 
#> Number of nonzero links: 230 
#> Percentage nonzero weights: 9.579342 
#> Average number of links: 4.693878 
#> 1 region with no links:
#> 21
#> 3 disjoint connected subgraphs
lw <- nb2listw(nb, style="B", zero.policy=TRUE)
joincount.test(HICRIME, lw)
#> 
#>  Join count test under nonfree sampling
#> 
#> data:  HICRIME 
#> weights: lw 
#> 
#> Std. deviate for low = 1.0036, p-value = 0.1578
#> alternative hypothesis: greater
#> sample estimates:
#> Same colour statistic           Expectation              Variance 
#>              35.00000              30.58511              19.35032 
#> 
#> 
#>  Join count test under nonfree sampling
#> 
#> data:  HICRIME 
#> weights: lw 
#> 
#> Std. deviate for high = 5.5716, p-value = 1.262e-08
#> alternative hypothesis: greater
#> sample estimates:
#> Same colour statistic           Expectation              Variance 
#>              52.00000              28.13830              18.34187 
#> 
```
