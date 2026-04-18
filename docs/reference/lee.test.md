# Lee's L test for spatial autocorrelation

Lee's L test for spatial autocorrelation using a spatial weights matrix
in weights list form. The assumptions underlying the test are sensitive
to the form of the graph of neighbour relationships and other factors,
and results may be checked against those of `lee.mc` permutations.

## Usage

``` r
lee.test(x, y, listw, zero.policy=attr(listw, "zero.policy"),
 alternative="greater", na.action=na.fail, spChk=NULL)
```

## Arguments

- x:

  a numeric vector the same length as the neighbours list in listw

- y:

  a numeric vector the same length as the neighbours list in listw

- listw:

  a `listw` object created for example by `nb2listw`

&nbsp;

- zero.policy:

  default `attr(listw, "zero.policy")` as set when `listw` was created,
  if attribute not set, use global option value; if TRUE assign zero to
  the lagged value of zones without neighbours, if FALSE assign NA

- alternative:

  a character string specifying the alternative hypothesis, must be one
  of greater (default), less or two.sided.

&nbsp;

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

## Value

A list with class `htest` containing the following components:

- statistic:

  the value of the standard deviate of Lee's L.

- p.value:

  the p-value of the test.

- estimate:

  the value of the observed Lee's L, its expectation and variance under
  the method assumption.

- alternative:

  a character string describing the alternative hypothesis.

- method:

  a character string giving the assumption used for calculating the
  standard deviate.

- data.name:

  a character string giving the name(s) of the data.

## Note

See Lee (2004) for details on how the asymptotic expectation and
variance of Lee's L is computed. In particular, check Lee (2004), table
1, page 1690.

This test may fail for large datasets as the computation of the
asymptotic expectation and variance requires the use of dense matrices.

## References

Lee (2004). A generalized significance testing method for global
measures of spatial association: an extension of the Mantel test.
Environment and Planning A 2004, volume 36, pages 1687 - 1703

## Author

Roger Bivand and Virgilio GÃ³mez-Rubio <Virgilio.Gomez@uclm.es>

## See also

[`lee`](https://r-spatial.github.io/spdep/reference/lee.md),
[`lee.mc`](https://r-spatial.github.io/spdep/reference/lee.mc.md),
[`listw2U`](https://r-spatial.github.io/spdep/reference/nb2listw.md)

## Examples

``` r
data(oldcol)
col.W <- nb2listw(COL.nb, style="W")
crime <- COL.OLD$CRIME

lee.test(crime, crime, col.W, zero.policy=TRUE)
#> 
#>  Lee's L statistic randomisation
#> 
#> data:  crime ,  crime 
#> weights: col.W  
#> 
#> Lee's L statistic standard deviate = 5.2343, p-value = 8.279e-08
#> alternative hypothesis: greater
#> sample estimates:
#> Lee's L statistic       Expectation          Variance 
#>       0.547064219       0.239417989       0.003454459 
#> 

#Test with missing values
x<-crime
y<-crime
x[1:5]<-NA
y[3:7]<-NA

lee.test(x, y, col.W, zero.policy=TRUE, na.action=na.omit)
#> 
#>  Lee's L statistic randomisation
#> 
#> data:  x ,  y 
#> weights: col.W 
#> omitted: 1, 2, 3, 4, 5, 6, 7 
#> 
#> Lee's L statistic standard deviate = 6.6873, p-value = 1.137e-11
#> alternative hypothesis: greater
#> sample estimates:
#> Lee's L statistic       Expectation          Variance 
#>       0.706469726       0.260143244       0.004454563 
#> 
#  lee.test(x, y, col.W, zero.policy=TRUE)#Stops with an error



data(boston, package="spData")
lw<-nb2listw(boston.soi)

x<-boston.c$CMEDV
y<-boston.c$CRIM

lee.test(x, y, lw, zero.policy=TRUE, alternative="less")
#> 
#>  Lee's L statistic randomisation
#> 
#> data:  x ,  y 
#> weights: lw  
#> 
#> Lee's L statistic standard deviate = -11.54, p-value < 2.2e-16
#> alternative hypothesis: less
#> sample estimates:
#> Lee's L statistic       Expectation          Variance 
#>      -0.326297206      -0.105040316       0.000367637 
#> 

#Test with missing values
x[1:5]<-NA
y[3:7]<-NA

lee.test(x, y, lw, zero.policy=TRUE, alternative="less", na.action=na.omit)
#> 
#>  Lee's L statistic randomisation
#> 
#> data:  x ,  y 
#> weights: lw 
#> omitted: 1, 2, 3, 4, 5, 6, 7 
#> 
#> Lee's L statistic standard deviate = -11.284, p-value < 2.2e-16
#> alternative hypothesis: less
#> sample estimates:
#> Lee's L statistic       Expectation          Variance 
#>     -0.3244748280     -0.1053154332      0.0003772109 
#> 













```
