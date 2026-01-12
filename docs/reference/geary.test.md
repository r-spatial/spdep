# Geary's C test for spatial autocorrelation

Geary's test for spatial autocorrelation using a spatial weights matrix
in weights list form. The assumptions underlying the test are sensitive
to the form of the graph of neighbour relationships and other factors,
and results may be checked against those of `geary.mc` permutations.

## Usage

``` r
geary.test(x, listw, randomisation=TRUE, zero.policy=attr(listw, "zero.policy"),
    alternative="greater", spChk=NULL, adjust.n=TRUE, na.action=na.fail,
    scale=TRUE)
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
  of "greater" (default), "less" or "two.sided".

- spChk:

  should the data vector names be checked against the spatial objects
  for identity integrity, TRUE, or FALSE, default NULL to use
  [`get.spChkOption()`](https://r-spatial.github.io/spdep/reference/set.spChkOption.md)

- adjust.n:

  default TRUE, if FALSE the number of observations is not adjusted for
  no-neighbour observations, if TRUE, the number of observations is
  adjusted

- na.action:

  a function (default `na.fail`), can also be `na.omit` or
  `na.exclude` - in these cases the weights list will be subsetted to
  remove NAs in the data. It may be necessary to set zero.policy to TRUE
  because this subsetting may create no-neighbour observations. Note
  that only weights lists created without using the glist argument to
  `nb2listw` may be subsetted. `na.pass` is not permitted.

- scale:

  default TRUE, may be FALSE to revert changes made to accommodate
  `localC` in November 2021 (see \#151)

## Value

A list with class `htest` containing the following components:

- statistic:

  the value of the standard deviate of Geary's C, in the order given in
  Cliff and Ord 1973, p. 21, which is (EC - C) / sqrt(VC), that is with
  the sign reversed with respect to the more usual (C - EC) / sqrt(VC);
  this means that the “greater” alternative for the Geary C test
  corresponds to the “greater” alternative for Moran's I test.

- p.value:

  the p-value of the test.

- estimate:

  the value of the observed Geary's C, its expectation and variance
  under the method assumption.

- alternative:

  a character string describing the alternative hypothesis.

- method:

  a character string giving the assumption used for calculating the
  standard deviate.

- data.name:

  a character string giving the name(s) of the data.

## References

Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 21, Cliff, A.
D., Ord, J. K. 1973 Spatial Autocorrelation, Pion, pp. 15-16, 21; Bivand
RS, Wong DWS 2018 Comparing implementations of global and local
indicators of spatial association. TEST, 27(3), 716–748
[doi:10.1007/s11749-018-0599-x](https://doi.org/10.1007/s11749-018-0599-x)

## Author

Roger Bivand <Roger.Bivand@nhh.no>

## Note

The derivation of the test (Cliff and Ord, 1981, p. 18) assumes that the
weights matrix is symmetric. For inherently non-symmetric matrices, such
as k-nearest neighbour matrices,
[`listw2U()`](https://r-spatial.github.io/spdep/reference/nb2listw.md)
can be used to make the matrix symmetric. In non-symmetric weights
matrix cases, the variance of the test statistic may be negative (thanks
to Franz Munoz I for a well documented bug report). Geary's C is
affected by non-symmetric weights under normality much more than Moran's
I. From 0.4-35, the sign of the standard deviate of C is changed to
match Cliff and Ord (1973, p. 21).

## See also

[`geary`](https://r-spatial.github.io/spdep/reference/geary.md),
[`geary.mc`](https://r-spatial.github.io/spdep/reference/geary.mc.md),
[`listw2U`](https://r-spatial.github.io/spdep/reference/nb2listw.md)

## Examples

``` r
data(oldcol)
geary.test(COL.OLD$CRIME, nb2listw(COL.nb, style="W"))
#> 
#>  Geary C test under randomisation
#> 
#> data:  COL.OLD$CRIME 
#> weights: nb2listw(COL.nb, style = "W")   
#> 
#> Geary C statistic standard deviate = 4.7605, p-value = 9.655e-07
#> alternative hypothesis: Expectation greater than statistic
#> sample estimates:
#> Geary C statistic       Expectation          Variance 
#>        0.52986993        1.00000000        0.00975278 
#> 
geary.test(COL.OLD$CRIME, nb2listw(COL.nb, style="W"),
 randomisation=FALSE)
#> 
#>  Geary C test under normality
#> 
#> data:  COL.OLD$CRIME 
#> weights: nb2listw(COL.nb, style = "W")   
#> 
#> Geary C statistic standard deviate = 4.6388, p-value = 1.752e-06
#> alternative hypothesis: Expectation greater than statistic
#> sample estimates:
#> Geary C statistic       Expectation          Variance 
#>        0.52986993        1.00000000        0.01027137 
#> 
colold.lags <- nblag(COL.nb, 3)
geary.test(COL.OLD$CRIME, nb2listw(colold.lags[[2]],
 style="W"))
#> 
#>  Geary C test under randomisation
#> 
#> data:  COL.OLD$CRIME 
#> weights: nb2listw(colold.lags[[2]], style = "W")   
#> 
#> Geary C statistic standard deviate = 2.2896, p-value = 0.01102
#> alternative hypothesis: Expectation greater than statistic
#> sample estimates:
#> Geary C statistic       Expectation          Variance 
#>       0.811285136       1.000000000       0.006793327 
#> 
geary.test(COL.OLD$CRIME, nb2listw(colold.lags[[3]],
 style="W"), alternative="greater")
#> 
#>  Geary C test under randomisation
#> 
#> data:  COL.OLD$CRIME 
#> weights: nb2listw(colold.lags[[3]], style = "W")   
#> 
#> Geary C statistic standard deviate = -1.5667, p-value = 0.9414
#> alternative hypothesis: Expectation greater than statistic
#> sample estimates:
#> Geary C statistic       Expectation          Variance 
#>       1.130277918       1.000000000       0.006914551 
#> 
print(is.symmetric.nb(COL.nb))
#> [1] TRUE
coords.OLD <- cbind(COL.OLD$X, COL.OLD$Y)
COL.k4.nb <- knn2nb(knearneigh(coords.OLD, 4))
print(is.symmetric.nb(COL.k4.nb))
#> [1] FALSE
geary.test(COL.OLD$CRIME, nb2listw(COL.k4.nb, style="W"))
#> 
#>  Geary C test under randomisation
#> 
#> data:  COL.OLD$CRIME 
#> weights: nb2listw(COL.k4.nb, style = "W")   
#> 
#> Geary C statistic standard deviate = 6.4415, p-value = 5.916e-11
#> alternative hypothesis: Expectation greater than statistic
#> sample estimates:
#> Geary C statistic       Expectation          Variance 
#>       0.399254423       1.000000000       0.008697812 
#> 
geary.test(COL.OLD$CRIME, nb2listw(COL.k4.nb, style="W"),
 randomisation=FALSE)
#> 
#>  Geary C test under normality
#> 
#> data:  COL.OLD$CRIME 
#> weights: nb2listw(COL.k4.nb, style = "W")   
#> 
#> Geary C statistic standard deviate = 6.2873, p-value = 1.615e-10
#> alternative hypothesis: Expectation greater than statistic
#> sample estimates:
#> Geary C statistic       Expectation          Variance 
#>       0.399254423       1.000000000       0.009129529 
#> 
cat("Note non-symmetric weights matrix - use listw2U()\n")
#> Note non-symmetric weights matrix - use listw2U()
geary.test(COL.OLD$CRIME, listw2U(nb2listw(COL.k4.nb,
 style="W")))
#> 
#>  Geary C test under randomisation
#> 
#> data:  COL.OLD$CRIME 
#> weights: listw2U(nb2listw(COL.k4.nb, style = "W"))   
#> 
#> Geary C statistic standard deviate = 6.4415, p-value = 5.916e-11
#> alternative hypothesis: Expectation greater than statistic
#> sample estimates:
#> Geary C statistic       Expectation          Variance 
#>       0.399254423       1.000000000       0.008697812 
#> 
geary.test(COL.OLD$CRIME, listw2U(nb2listw(COL.k4.nb,
 style="W")), randomisation=FALSE)
#> 
#>  Geary C test under normality
#> 
#> data:  COL.OLD$CRIME 
#> weights: listw2U(nb2listw(COL.k4.nb, style = "W"))   
#> 
#> Geary C statistic standard deviate = 6.2873, p-value = 1.615e-10
#> alternative hypothesis: Expectation greater than statistic
#> sample estimates:
#> Geary C statistic       Expectation          Variance 
#>       0.399254423       1.000000000       0.009129529 
#> 
crime <- COL.OLD$CRIME
is.na(crime) <- sample(1:length(crime), 10)
try(geary.test(crime, nb2listw(COL.nb, style="W"), na.action=na.fail))
#> Error in na.fail.default(x) : missing values in object
geary.test(crime, nb2listw(COL.nb, style="W"), zero.policy=TRUE,
 na.action=na.omit)
#> 
#>  Geary C test under randomisation
#> 
#> data:  crime 
#> weights: nb2listw(COL.nb, style = "W") 
#> omitted: 3, 4, 10, 13, 21, 23, 27, 29, 31, 38  
#> 
#> Geary C statistic standard deviate = 4.1449, p-value = 1.7e-05
#> alternative hypothesis: Expectation greater than statistic
#> sample estimates:
#> Geary C statistic       Expectation          Variance 
#>        0.46652084        1.00000000        0.01656579 
#> 
geary.test(crime, nb2listw(COL.nb, style="W"), zero.policy=TRUE,
 na.action=na.exclude)
#> 
#>  Geary C test under randomisation
#> 
#> data:  crime 
#> weights: nb2listw(COL.nb, style = "W") 
#> omitted: 3, 4, 10, 13, 21, 23, 27, 29, 31, 38  
#> 
#> Geary C statistic standard deviate = 4.1449, p-value = 1.7e-05
#> alternative hypothesis: Expectation greater than statistic
#> sample estimates:
#> Geary C statistic       Expectation          Variance 
#>        0.46652084        1.00000000        0.01656579 
#> 
try(geary.test(crime, nb2listw(COL.nb, style="W"), na.action=na.pass))
#> Error in geary.test(crime, nb2listw(COL.nb, style = "W"), na.action = na.pass) : 
#>   na.pass not permitted
```
