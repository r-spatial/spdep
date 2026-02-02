# Local Moran's I statistic

The local spatial statistic Moran's I is calculated for each zone based
on the spatial weights object used. The values returned include a
Z-value, and may be used as a diagnostic tool. The statistic and its
expectation and variance were given in Anselin (1995), but those from
Sokal et al. (1998) are implemented here.

## Usage

``` r
localmoran(x, listw, zero.policy=attr(listw, "zero.policy"), na.action=na.fail,
        conditional=TRUE, alternative = "two.sided", mlvar=TRUE,
        spChk=NULL, adjust.x=FALSE)
localmoran_perm(x, listw, nsim=499, zero.policy=attr(listw, "zero.policy"), 
        na.action=na.fail, alternative = "two.sided", mlvar=TRUE,
        spChk=NULL, adjust.x=FALSE, sample_Ei=TRUE, iseed=NULL,
        no_repeat_in_row=FALSE)
```

## Arguments

- x:

  a numeric vector the same length as the neighbours list in listw

- listw:

  a `listw` object created for example by `nb2listw`

- zero.policy:

  default default `attr(listw, "zero.policy")` as set when `listw` was
  created, if attribute not set, use global option value; if TRUE assign
  zero to the lagged value of zones without neighbours, if FALSE assign
  NA

- na.action:

  a function (default `na.fail`), can also be `na.omit` or
  `na.exclude` - in these cases the weights list will be subsetted to
  remove NAs in the data. It may be necessary to set zero.policy to TRUE
  because this subsetting may create no-neighbour observations. Note
  that only weights lists created without using the glist argument to
  `nb2listw` may be subsetted. If `na.pass` is used, zero is substituted
  for NA values in calculating the spatial lag. (Note that na.exclude
  will only work properly starting from R 1.9.0, na.omit and na.exclude
  assign the wrong classes in 1.8.\*)

- conditional:

  default TRUE: expectation and variance are calculated using the
  conditional randomization null (Sokal 1998 Eqs. A7 & A8). Elaboration
  of these changes available in Sauer et al. (2021). If FALSE:
  expectation and variance are calculated using the total randomization
  null (Sokal 1998 Eqs. A3 & A4).

- alternative:

  a character string specifying the alternative hypothesis, must be one
  of greater, less or two.sided (default).

- mlvar:

  default TRUE: values of local Moran's I are reported using the
  variance of the variable of interest (sum of squared deviances over
  n), but can be reported as the sample variance, dividing by (n-1)
  instead; both are used in other implementations.

- spChk:

  should the data vector names be checked against the spatial objects
  for identity integrity, TRUE, or FALSE, default NULL to use
  [`get.spChkOption()`](https://r-spatial.github.io/spdep/reference/set.spChkOption.md)

- adjust.x:

  default FALSE, if TRUE, x values of observations with no neighbours
  are omitted in the mean of x

- nsim:

  default 499, number of conditonal permutation simulations

- sample_Ei:

  default TRUE; if conditional permutation, use the sample \$E_i\$
  values, or the analytical values, leaving only variances calculated by
  simulation.

- iseed:

  default NULL, used to set the seed; the output will only be
  reproducible if the count of CPU cores across which computation is
  distributed is the same

- no_repeat_in_row:

  default `FALSE`, if `TRUE`, sample conditionally in each row without
  replacements to avoid duplicate values,
  <https://github.com/r-spatial/spdep/issues/124>

## Details

The values of local Moran's I are divided by the variance (or sample
variance) of the variable of interest to accord with Table 1, p. 103,
and formula (12), p. 99, in Anselin (1995), rather than his formula (7),
p. 98. The variance of the local Moran statistic is taken from Sokal et
al. (1998) p. 334, equations 4 & 5 or equations 7 & 8 located depending
on user specification. By default, the implementation divides by n, not
(n-1) in calculating the variance and higher moments. Conditional code
contributed by Jeff Sauer and Levi Wolf.

## Value

- Ii:

  local moran statistic

- E.Ii:

  expectation of local moran statistic; for `localmoran_perm`the
  permutation sample means

- Var.Ii:

  variance of local moran statistic; for `localmoran_perm`the
  permutation sample standard deviations

- Z.Ii:

  standard deviate of local moran statistic; for `localmoran_perm` based
  on permutation sample means and standard deviations

- Pr():

  p-value of local moran statistic using
  [`pnorm()`](https://rdrr.io/r/stats/Normal.html); for
  `localmoran_perm` using standard deviatse based on permutation sample
  means and standard deviations

- Pr() Sim:

  For `localmoran_perm`, [`rank()`](https://rdrr.io/r/base/rank.html)
  and [`punif()`](https://rdrr.io/r/stats/Uniform.html) of observed
  statistic rank for \[0, 1\] p-values using `alternative=`

- Pr(folded) Sim:

  the simulation folded \[0, 0.5\] range ranked p-value (based on
  <https://github.com/pysal/esda/blob/4a63e0b5df1e754b17b5f1205b8cadcbecc5e061/esda/crand.py#L211-L213>)

- Skewness:

  For `localmoran_perm`, the output of
  [`e1071::skewness()`](https://rdrr.io/pkg/e1071/man/skewness.html) for
  the permutation samples underlying the standard deviates

- Kurtosis:

  For `localmoran_perm`, the output of
  [`e1071::kurtosis()`](https://rdrr.io/pkg/e1071/man/kurtosis.html) for
  the permutation samples underlying the standard deviates

In addition, an attribute data frame `"quadr"` with mean and median
quadrant columns, and a column splitting on the demeaned variable and
lagged demeaned variable at zero.

## Note

Conditional permutations added for comparative purposes; permutations
are over the whole data vector omitting the observation itself. For
p-value adjustment, use
[`p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) or
[`p.adjustSP()`](https://r-spatial.github.io/spdep/reference/p.adjustSP.md)
on the output vector.

## References

Anselin, L. 1995. Local indicators of spatial association, Geographical
Analysis, 27, 93–115; Getis, A. and Ord, J. K. 1996 Local spatial
statistics: an overview. In P. Longley and M. Batty (eds) *Spatial
analysis: modelling in a GIS environment* (Cambridge: Geoinformation
International), 261–277; Sokal, R. R, Oden, N. L. and Thomson, B. A.
1998. Local Spatial Autocorrelation in a Biological Model. Geographical
Analysis, 30. 331–354; Bivand RS, Wong DWS 2018 Comparing
implementations of global and local indicators of spatial association.
TEST, 27(3), 716–748
[doi:10.1007/s11749-018-0599-x](https://doi.org/10.1007/s11749-018-0599-x)
; Sauer, J., Oshan, T. M., Rey, S., & Wolf, L. J. 2021. The Importance
of Null Hypotheses: Understanding Differences in Local Moran’s under
Heteroskedasticity. Geographical Analysis.
[doi:10.1111/gean.12304](https://doi.org/10.1111/gean.12304)

Bivand, R. (2022), R Packages for Analyzing Spatial Data: A Comparative
Case Study with Areal Data. Geographical Analysis, 54(3), 488-518.
[doi:10.1111/gean.12319](https://doi.org/10.1111/gean.12319)

## Author

Roger Bivand <Roger.Bivand@nhh.no>

## See also

[`localG`](https://r-spatial.github.io/spdep/reference/localG.md)

## Examples

``` r
data(afcon, package="spData")
oid <- order(afcon$id)
resI <- localmoran(afcon$totcon, nb2listw(paper.nb))
printCoefmat(data.frame(resI[oid,], row.names=afcon$name[oid]),
 check.names=FALSE)
#>                                   Ii        E.Ii      Var.Ii        Z.Ii
#> THE GAMBIA                3.7523e-01 -2.4317e-02  9.9646e-01  4.0025e-01
#> MALI                      4.6363e-01 -2.1841e-02  1.0896e-01  1.4708e+00
#> SENEGAL                   2.5670e-01 -3.4444e-03  3.3339e-02  1.4248e+00
#> BENIN                     1.9412e-01 -2.4556e-03  2.3792e-02  1.2744e+00
#> MAURITANIA                9.7053e-02 -5.7508e-03  5.5533e-02  4.3625e-01
#> NIGER                     2.3071e-01 -1.9459e-02  9.7310e-02  8.0198e-01
#> IVORY COAST               2.9004e-01 -6.9359e-03  5.2072e-02  1.3014e+00
#> GUINEA                    1.8263e-01 -2.2246e-03  1.6780e-02  1.4270e+00
#> BURKINA FASO              5.0828e-01 -1.9893e-02  1.1942e-01  1.5284e+00
#> LIBERIA                   1.8565e-01 -2.7127e-03  3.5982e-02  9.9300e-01
#> SIERRA LEONE              2.6523e-01 -1.6994e-02  3.4204e-01  4.8257e-01
#> GHANA                     1.4764e-01 -1.3414e-03  1.7817e-02  1.1161e+00
#> TOGO                      2.1934e-01 -4.9892e-03  6.6025e-02  8.7305e-01
#> CAMEROON                  2.5925e-01 -1.1009e-02  8.2313e-02  9.4198e-01
#> NIGERIA                   1.1377e-01 -9.6126e-04  9.3272e-03  1.1880e+00
#> GABON                     2.0366e-01 -5.4771e-03  1.1153e-01  6.2625e-01
#> CENTRAL AFRICAN REPUBLIC -4.4206e-01 -1.0600e-02  7.9287e-02 -1.5323e+00
#> CHAD                     -1.0528e-01 -4.0998e-03  2.5008e-02 -6.3985e-01
#> CONGO                     1.1380e-02 -8.5953e-04  8.3410e-03  1.3402e-01
#> ZAIRE                     7.0978e-01 -5.9545e-02  2.0906e-01  1.6826e+00
#> ANGOLA                    1.1797e-01 -6.2140e-04  8.2594e-03  1.3050e+00
#> UGANDA                    1.9425e+00 -6.2812e-02  4.4503e-01  3.0060e+00
#> KENYA                     1.1969e+00 -1.6803e-02  1.2489e-01  3.4344e+00
#> TANZANIA                  2.7185e-01 -4.6254e-02  1.9107e-01  7.2774e-01
#> BURUNDI                  -4.8428e-01 -1.1009e-02  1.4481e-01 -1.2437e+00
#> RWANDA                   -7.5236e-01 -1.4730e-02  1.4096e-01 -1.9647e+00
#> SOMALIA                   4.5277e-01 -1.1751e-02  2.3778e-01  9.5260e-01
#> ETHIOPIA                  7.2512e-01 -5.4929e-03  7.2655e-02  2.7106e+00
#> ZAMBIA                    4.2160e-02 -8.1691e-04  3.5354e-03  7.2280e-01
#> ZIMBABWE                 -9.5068e-03 -6.0969e-03  5.8855e-02 -1.4056e-02
#> MALAWI                   -2.2888e-01 -1.0284e-02  1.3537e-01 -5.9413e-01
#> MOZAMBIQUE                1.6790e-02 -6.1629e-03  3.7515e-02  1.1850e-01
#> SOUTH AFRICA             -1.8254e-01 -5.4306e-03  2.7546e-02 -1.0671e+00
#> LESOTHO                  -4.1935e-01 -1.9263e-02  7.9348e-01 -4.4914e-01
#> BOTSWANA                 -3.9316e-03 -1.4141e-04  1.8805e-03 -8.7403e-02
#> SWAZILAND                 1.6684e-02 -2.8611e-02  5.6905e-01  6.0045e-02
#> MOROCCO                  -9.6961e-02 -5.1445e-03  1.0479e-01 -2.8363e-01
#> ALGERIA                  -1.0037e-02 -9.7828e-05  5.9914e-04 -4.0605e-01
#> TUNISIA                   5.3873e-03 -3.0273e-06  6.1985e-05  6.8466e-01
#> LIBYA                     8.0382e-01 -1.9923e-02  1.1960e-01  2.3820e+00
#> SUDAN                     2.9878e+00 -2.2835e-01  7.6320e-01  3.6814e+00
#> EGYPT                     6.9467e+00 -2.9968e-01  4.2971e+00  3.4957e+00
#>                          Pr.z....E.Ii..
#> THE GAMBIA                       0.6890
#> MALI                             0.1414
#> SENEGAL                          0.1542
#> BENIN                            0.2025
#> MAURITANIA                       0.6627
#> NIGER                            0.4226
#> IVORY COAST                      0.1931
#> GUINEA                           0.1536
#> BURKINA FASO                     0.1264
#> LIBERIA                          0.3207
#> SIERRA LEONE                     0.6294
#> GHANA                            0.2644
#> TOGO                             0.3826
#> CAMEROON                         0.3462
#> NIGERIA                          0.2348
#> GABON                            0.5312
#> CENTRAL AFRICAN REPUBLIC         0.1255
#> CHAD                             0.5223
#> CONGO                            0.8934
#> ZAIRE                            0.0925
#> ANGOLA                           0.1919
#> UGANDA                           0.0026
#> KENYA                            0.0006
#> TANZANIA                         0.4668
#> BURUNDI                          0.2136
#> RWANDA                           0.0494
#> SOMALIA                          0.3408
#> ETHIOPIA                         0.0067
#> ZAMBIA                           0.4698
#> ZIMBABWE                         0.9888
#> MALAWI                           0.5524
#> MOZAMBIQUE                       0.9057
#> SOUTH AFRICA                     0.2859
#> LESOTHO                          0.6533
#> BOTSWANA                         0.9304
#> SWAZILAND                        0.9521
#> MOROCCO                          0.7767
#> ALGERIA                          0.6847
#> TUNISIA                          0.4936
#> LIBYA                            0.0172
#> SUDAN                            0.0002
#> EGYPT                            0.0005
hist(resI[,5])

mean(resI[,1])
#> [1] 0.4167956
sum(resI[,1])/Szero(nb2listw(paper.nb))
#> [1] 0.4167956
moran.test(afcon$totcon, nb2listw(paper.nb))
#> 
#>  Moran I test under randomisation
#> 
#> data:  afcon$totcon  
#> weights: nb2listw(paper.nb)    
#> 
#> Moran I statistic standard deviate = 4.3485, p-value = 6.854e-06
#> alternative hypothesis: greater
#> sample estimates:
#> Moran I statistic       Expectation          Variance 
#>        0.41679563       -0.02439024        0.01029358 
#> 
# note equality for mean() only when the sum of weights equals
# the number of observations (thanks to Juergen Symanzik)
resI <- localmoran(afcon$totcon, nb2listw(paper.nb))
printCoefmat(data.frame(resI[oid,], row.names=afcon$name[oid]),
 check.names=FALSE)
#>                                   Ii        E.Ii      Var.Ii        Z.Ii
#> THE GAMBIA                3.7523e-01 -2.4317e-02  9.9646e-01  4.0025e-01
#> MALI                      4.6363e-01 -2.1841e-02  1.0896e-01  1.4708e+00
#> SENEGAL                   2.5670e-01 -3.4444e-03  3.3339e-02  1.4248e+00
#> BENIN                     1.9412e-01 -2.4556e-03  2.3792e-02  1.2744e+00
#> MAURITANIA                9.7053e-02 -5.7508e-03  5.5533e-02  4.3625e-01
#> NIGER                     2.3071e-01 -1.9459e-02  9.7310e-02  8.0198e-01
#> IVORY COAST               2.9004e-01 -6.9359e-03  5.2072e-02  1.3014e+00
#> GUINEA                    1.8263e-01 -2.2246e-03  1.6780e-02  1.4270e+00
#> BURKINA FASO              5.0828e-01 -1.9893e-02  1.1942e-01  1.5284e+00
#> LIBERIA                   1.8565e-01 -2.7127e-03  3.5982e-02  9.9300e-01
#> SIERRA LEONE              2.6523e-01 -1.6994e-02  3.4204e-01  4.8257e-01
#> GHANA                     1.4764e-01 -1.3414e-03  1.7817e-02  1.1161e+00
#> TOGO                      2.1934e-01 -4.9892e-03  6.6025e-02  8.7305e-01
#> CAMEROON                  2.5925e-01 -1.1009e-02  8.2313e-02  9.4198e-01
#> NIGERIA                   1.1377e-01 -9.6126e-04  9.3272e-03  1.1880e+00
#> GABON                     2.0366e-01 -5.4771e-03  1.1153e-01  6.2625e-01
#> CENTRAL AFRICAN REPUBLIC -4.4206e-01 -1.0600e-02  7.9287e-02 -1.5323e+00
#> CHAD                     -1.0528e-01 -4.0998e-03  2.5008e-02 -6.3985e-01
#> CONGO                     1.1380e-02 -8.5953e-04  8.3410e-03  1.3402e-01
#> ZAIRE                     7.0978e-01 -5.9545e-02  2.0906e-01  1.6826e+00
#> ANGOLA                    1.1797e-01 -6.2140e-04  8.2594e-03  1.3050e+00
#> UGANDA                    1.9425e+00 -6.2812e-02  4.4503e-01  3.0060e+00
#> KENYA                     1.1969e+00 -1.6803e-02  1.2489e-01  3.4344e+00
#> TANZANIA                  2.7185e-01 -4.6254e-02  1.9107e-01  7.2774e-01
#> BURUNDI                  -4.8428e-01 -1.1009e-02  1.4481e-01 -1.2437e+00
#> RWANDA                   -7.5236e-01 -1.4730e-02  1.4096e-01 -1.9647e+00
#> SOMALIA                   4.5277e-01 -1.1751e-02  2.3778e-01  9.5260e-01
#> ETHIOPIA                  7.2512e-01 -5.4929e-03  7.2655e-02  2.7106e+00
#> ZAMBIA                    4.2160e-02 -8.1691e-04  3.5354e-03  7.2280e-01
#> ZIMBABWE                 -9.5068e-03 -6.0969e-03  5.8855e-02 -1.4056e-02
#> MALAWI                   -2.2888e-01 -1.0284e-02  1.3537e-01 -5.9413e-01
#> MOZAMBIQUE                1.6790e-02 -6.1629e-03  3.7515e-02  1.1850e-01
#> SOUTH AFRICA             -1.8254e-01 -5.4306e-03  2.7546e-02 -1.0671e+00
#> LESOTHO                  -4.1935e-01 -1.9263e-02  7.9348e-01 -4.4914e-01
#> BOTSWANA                 -3.9316e-03 -1.4141e-04  1.8805e-03 -8.7403e-02
#> SWAZILAND                 1.6684e-02 -2.8611e-02  5.6905e-01  6.0045e-02
#> MOROCCO                  -9.6961e-02 -5.1445e-03  1.0479e-01 -2.8363e-01
#> ALGERIA                  -1.0037e-02 -9.7828e-05  5.9914e-04 -4.0605e-01
#> TUNISIA                   5.3873e-03 -3.0273e-06  6.1985e-05  6.8466e-01
#> LIBYA                     8.0382e-01 -1.9923e-02  1.1960e-01  2.3820e+00
#> SUDAN                     2.9878e+00 -2.2835e-01  7.6320e-01  3.6814e+00
#> EGYPT                     6.9467e+00 -2.9968e-01  4.2971e+00  3.4957e+00
#>                          Pr.z....E.Ii..
#> THE GAMBIA                       0.6890
#> MALI                             0.1414
#> SENEGAL                          0.1542
#> BENIN                            0.2025
#> MAURITANIA                       0.6627
#> NIGER                            0.4226
#> IVORY COAST                      0.1931
#> GUINEA                           0.1536
#> BURKINA FASO                     0.1264
#> LIBERIA                          0.3207
#> SIERRA LEONE                     0.6294
#> GHANA                            0.2644
#> TOGO                             0.3826
#> CAMEROON                         0.3462
#> NIGERIA                          0.2348
#> GABON                            0.5312
#> CENTRAL AFRICAN REPUBLIC         0.1255
#> CHAD                             0.5223
#> CONGO                            0.8934
#> ZAIRE                            0.0925
#> ANGOLA                           0.1919
#> UGANDA                           0.0026
#> KENYA                            0.0006
#> TANZANIA                         0.4668
#> BURUNDI                          0.2136
#> RWANDA                           0.0494
#> SOMALIA                          0.3408
#> ETHIOPIA                         0.0067
#> ZAMBIA                           0.4698
#> ZIMBABWE                         0.9888
#> MALAWI                           0.5524
#> MOZAMBIQUE                       0.9057
#> SOUTH AFRICA                     0.2859
#> LESOTHO                          0.6533
#> BOTSWANA                         0.9304
#> SWAZILAND                        0.9521
#> MOROCCO                          0.7767
#> ALGERIA                          0.6847
#> TUNISIA                          0.4936
#> LIBYA                            0.0172
#> SUDAN                            0.0002
#> EGYPT                            0.0005
hist(p.adjust(resI[,5], method="bonferroni"))

totcon <-afcon$totcon
is.na(totcon) <- sample(1:length(totcon), 5)
totcon
#>  [1] 1363 1421   NA   NA 5246  811  299  358  895 4751 1878  933  347 1130  241
#> [16]  604 1015  998 2122 1090  848  618   NA  423  980 3087 2273 3134 1142  824
#> [31] 2881   NA  604 1528 1554   NA  792  795 1266 1875  147  363
resI.na <- localmoran(totcon, nb2listw(paper.nb), na.action=na.exclude,
 zero.policy=TRUE)
if (class(attr(resI.na, "na.action")) == "exclude") {
 print(data.frame(resI.na[oid,], row.names=afcon$name[oid]), digits=2)
} else print(resI.na, digits=2)
#>                                Ii     E.Ii  Var.Ii   Z.Ii Pr.z....E.Ii..
#> THE GAMBIA                0.37105 -2.7e-02 9.6e-01  0.406        0.68494
#> MALI                      0.44799 -2.4e-02 1.2e-01  1.341        0.17986
#> SENEGAL                   0.25571 -4.0e-03 3.4e-02  1.418        0.15618
#> BENIN                     0.19536 -2.9e-03 2.4e-02  1.272        0.20348
#> MAURITANIA                0.20397 -6.5e-03 7.5e-02  0.767        0.44324
#> NIGER                     0.39719 -2.1e-02 1.1e-01  1.257        0.20865
#> IVORY COAST                    NA       NA      NA     NA             NA
#> GUINEA                    0.18960 -2.6e-03 2.2e-02  1.292        0.19628
#> BURKINA FASO              0.50165 -2.2e-02 1.4e-01  1.397        0.16235
#> LIBERIA                   0.19068 -3.2e-03 5.7e-02  0.814        0.41584
#> SIERRA LEONE              0.26508 -1.9e-02 3.3e-01  0.493        0.62170
#> GHANA                     0.16227 -1.6e-03 2.9e-02  0.959        0.33762
#> TOGO                      0.21902 -5.7e-03 6.6e-02  0.876        0.38087
#> CAMEROON                  0.25806 -1.2e-02 7.9e-02  0.959        0.33744
#> NIGERIA                   0.11801 -1.2e-03 1.0e-02  1.187        0.23527
#> GABON                     0.20388 -6.2e-03 1.1e-01  0.630        0.52856
#> CENTRAL AFRICAN REPUBLIC -0.41241 -1.2e-02 7.7e-02 -1.448        0.14772
#> CHAD                     -0.04424 -4.7e-03 3.1e-02 -0.226        0.82148
#> CONGO                     0.01460 -1.1e-03 9.1e-03  0.164        0.86947
#> ZAIRE                     0.85429 -6.2e-02 2.2e-01  1.976        0.04820
#> ANGOLA                    0.09680 -5.3e-04 6.2e-03  1.236        0.21649
#> UGANDA                    2.50702 -6.5e-02 5.2e-01  3.578        0.00035
#> KENYA                     1.08291 -1.7e-02 1.1e-01  3.308        0.00094
#> TANZANIA                  0.61571 -4.8e-02 2.4e-01  1.350        0.17690
#> BURUNDI                  -0.93318 -1.2e-02 2.2e-01 -1.973        0.04845
#> RWANDA                         NA       NA      NA     NA             NA
#> SOMALIA                   0.40246 -1.2e-02 2.1e-01  0.901        0.36766
#> ETHIOPIA                  0.64672 -5.4e-03 6.3e-02  2.598        0.00937
#> ZAMBIA                    0.05293 -7.2e-04 3.2e-03  0.955        0.33980
#> ZIMBABWE                 -0.00139 -6.9e-03 5.8e-02  0.023        0.98173
#> MALAWI                         NA       NA      NA     NA             NA
#> MOZAMBIQUE               -0.03566 -7.0e-03 4.5e-02 -0.135        0.89292
#> SOUTH AFRICA             -0.17136 -5.4e-03 2.3e-02 -1.084        0.27841
#> LESOTHO                  -0.38478 -2.1e-02 7.7e-01 -0.415        0.67835
#> BOTSWANA                 -0.00306 -2.2e-04 2.6e-03 -0.056        0.95534
#> SWAZILAND                 0.03234 -3.1e-02 5.5e-01  0.086        0.93130
#> MOROCCO                        NA       NA      NA     NA             NA
#> ALGERIA                  -0.02619 -5.8e-05 4.9e-04 -1.182        0.23702
#> TUNISIA                  -0.00022 -6.4e-07 2.4e-05 -0.045        0.96423
#> LIBYA                          NA       NA      NA     NA             NA
#> SUDAN                     2.75689 -2.4e-01 8.0e-01  3.353        0.00080
#> EGYPT                     9.90940 -3.2e-01 8.0e+00  3.617        0.00030
resG <- localG(afcon$totcon, nb2listw(include.self(paper.nb)))
print(data.frame(resG[oid], row.names=afcon$name[oid]), digits=2)
#>                          resG.oid.
#> THE GAMBIA                  -0.984
#> MALI                        -1.699
#> SENEGAL                     -1.463
#> BENIN                       -1.301
#> MAURITANIA                  -0.605
#> NIGER                       -1.049
#> IVORY COAST                 -1.417
#> GUINEA                      -1.449
#> BURKINA FASO                -1.751
#> LIBERIA                     -1.041
#> SIERRA LEONE                -0.870
#> GHANA                       -1.103
#> TOGO                        -0.991
#> CAMEROON                    -1.133
#> NIGERIA                     -1.173
#> GABON                       -0.789
#> CENTRAL AFRICAN REPUBLIC     1.173
#> CHAD                         0.463
#> CONGO                       -0.203
#> ZAIRE                        2.023
#> ANGOLA                       1.235
#> UGANDA                       3.336
#> KENYA                        3.503
#> TANZANIA                     1.098
#> BURUNDI                      0.774
#> RWANDA                       1.457
#> SOMALIA                      1.183
#> ETHIOPIA                     2.627
#> ZAMBIA                       0.753
#> ZIMBABWE                    -0.200
#> MALAWI                       0.212
#> MOZAMBIQUE                  -0.288
#> SOUTH AFRICA                -0.868
#> LESOTHO                     -0.298
#> BOTSWANA                     0.041
#> SWAZILAND                   -0.659
#> MOROCCO                      0.022
#> ALGERIA                     -0.363
#> TUNISIA                      0.579
#> LIBYA                        2.553
#> SUDAN                        4.039
#> EGYPT                        4.421
set.seed(1)
resI_p <- localmoran_perm(afcon$totcon, nb2listw(paper.nb))
printCoefmat(data.frame(resI_p[oid,], row.names=afcon$name[oid]),
 check.names=FALSE)
#>                                   Ii        E.Ii      Var.Ii        Z.Ii
#> THE GAMBIA                3.7523e-01 -3.0853e-02  9.0476e-01  4.2692e-01
#> MALI                      4.6363e-01 -2.1690e-02  1.2192e-01  1.3899e+00
#> SENEGAL                   2.5670e-01 -1.1288e-02  3.2719e-02  1.4816e+00
#> BENIN                     1.9412e-01 -6.5012e-03  2.6046e-02  1.2431e+00
#> MAURITANIA                9.7053e-02 -2.0574e-03  5.7209e-02  4.1437e-01
#> NIGER                     2.3071e-01 -2.8428e-02  1.1981e-01  7.4867e-01
#> IVORY COAST               2.9004e-01 -8.1011e-03  5.5564e-02  1.2648e+00
#> GUINEA                    1.8263e-01  5.3454e-04  1.8397e-02  1.3425e+00
#> BURKINA FASO              5.0828e-01 -2.9727e-02  1.4019e-01  1.4369e+00
#> LIBERIA                   1.8565e-01 -1.9247e-03  3.6524e-02  9.8147e-01
#> SIERRA LEONE              2.6523e-01  1.0225e-02  3.1505e-01  4.5432e-01
#> GHANA                     1.4764e-01 -6.7639e-03  1.8646e-02  1.1307e+00
#> TOGO                      2.1934e-01 -1.8012e-03  6.5187e-02  8.6616e-01
#> CAMEROON                  2.5925e-01 -1.0726e-02  9.4789e-02  8.7688e-01
#> NIGERIA                   1.1377e-01 -2.5672e-03  8.7720e-03  1.2421e+00
#> GABON                     2.0366e-01 -3.5024e-03  1.0958e-01  6.2584e-01
#> CENTRAL AFRICAN REPUBLIC -4.4206e-01  2.3901e-02  8.1700e-02 -1.6302e+00
#> CHAD                     -1.0528e-01  5.9042e-03  3.0168e-02 -6.4016e-01
#> CONGO                     1.1380e-02  7.0349e-03  8.9435e-03  4.5949e-02
#> ZAIRE                     7.0978e-01 -2.9589e-02  2.9068e-01  1.3714e+00
#> ANGOLA                    1.1797e-01 -3.0421e-04  8.6526e-03  1.2716e+00
#> UGANDA                    1.9425e+00 -7.0491e-02  5.2046e-01  2.7903e+00
#> KENYA                     1.1969e+00 -1.6362e-02  1.3644e-01  3.2847e+00
#> TANZANIA                  2.7185e-01 -7.3309e-02  2.3454e-01  7.1272e-01
#> BURUNDI                  -4.8428e-01  9.5351e-03  1.4866e-01 -1.2807e+00
#> RWANDA                   -7.5236e-01 -4.0426e-02  1.7628e-01 -1.6956e+00
#> SOMALIA                   4.5277e-01 -1.1165e-02  2.6266e-01  9.0523e-01
#> ETHIOPIA                  7.2512e-01  8.3618e-04  7.8477e-02  2.5855e+00
#> ZAMBIA                    4.2160e-02 -5.1840e-03  3.6852e-03  7.7989e-01
#> ZIMBABWE                 -9.5068e-03 -6.5338e-03  6.5357e-02 -1.1629e-02
#> MALAWI                   -2.2888e-01  3.6276e-03  1.3417e-01 -6.3477e-01
#> MOZAMBIQUE                1.6790e-02 -2.0694e-03  4.5436e-02  8.8476e-02
#> SOUTH AFRICA             -1.8254e-01 -1.1611e-02  3.0240e-02 -9.8293e-01
#> LESOTHO                  -4.1935e-01 -6.6258e-02  8.7314e-01 -3.7787e-01
#> BOTSWANA                 -3.9316e-03 -4.7713e-03  2.2167e-03  1.7833e-02
#> SWAZILAND                 1.6684e-02 -8.4486e-02  6.4807e-01  1.2567e-01
#> MOROCCO                  -9.6961e-02 -2.1623e-03  1.1253e-01 -2.8260e-01
#> ALGERIA                  -1.0037e-02 -3.4326e-04  6.3541e-04 -3.8456e-01
#> TUNISIA                   5.3873e-03 -1.6456e-04  6.4187e-05  6.9298e-01
#> LIBYA                     8.0382e-01 -2.9038e-02  1.3319e-01  2.2821e+00
#> SUDAN                     2.9878e+00 -2.5952e-01  8.2522e-01  3.5747e+00
#> EGYPT                     6.9467e+00 -3.7869e-01  4.1318e+00  3.6038e+00
#>                          Pr.z....E.Ii.. Pr.z....E.Ii...Sim Pr.folded..Sim
#> THE GAMBIA                   6.6944e-01         8.3200e-01     4.2800e-01
#> MALI                         1.6455e-01         8.8000e-02     4.4000e-02
#> SENEGAL                      1.3846e-01         6.4000e-02     3.2000e-02
#> BENIN                        2.1384e-01         1.3200e-01     6.6000e-02
#> MAURITANIA                   6.7860e-01         7.5600e-01     3.7800e-01
#> NIGER                        4.5406e-01         5.0000e-01     2.5000e-01
#> IVORY COAST                  2.0595e-01         1.1600e-01     5.8000e-02
#> GUINEA                       1.7943e-01         7.2000e-02     3.6000e-02
#> BURKINA FASO                 1.5075e-01         9.2000e-02     4.6000e-02
#> LIBERIA                      3.2636e-01         3.0400e-01     1.5200e-01
#> SIERRA LEONE                 6.4960e-01         8.2400e-01     4.1200e-01
#> GHANA                        2.5816e-01         1.8400e-01     9.2000e-02
#> TOGO                         3.8640e-01         3.6800e-01     1.8400e-01
#> CAMEROON                     3.8055e-01         4.3200e-01     2.1600e-01
#> NIGERIA                      2.1418e-01         1.4400e-01     7.2000e-02
#> GABON                        5.3142e-01         5.8400e-01     2.9600e-01
#> CENTRAL AFRICAN REPUBLIC     1.0306e-01         1.8400e-01     9.2000e-02
#> CHAD                         5.2207e-01         4.5200e-01     2.2600e-01
#> CONGO                        9.6335e-01         8.5600e-01     4.2800e-01
#> ZAIRE                        1.7026e-01         2.0000e-01     1.0000e-01
#> ANGOLA                       2.0353e-01         2.3600e-01     1.1800e-01
#> UGANDA                       5.2662e-03         4.0000e-02     2.0000e-02
#> KENYA                        1.0210e-03         8.0000e-03     4.0000e-03
#> TANZANIA                     4.7602e-01         4.2400e-01     2.1200e-01
#> BURUNDI                      2.0029e-01         2.6800e-01     1.3400e-01
#> RWANDA                       8.9953e-02         1.4000e-01     7.0000e-02
#> SOMALIA                      3.6535e-01         2.8800e-01     1.4600e-01
#> ETHIOPIA                     9.7244e-03         3.6000e-02     1.8000e-02
#> ZAMBIA                       4.3546e-01         4.1200e-01     2.0600e-01
#> ZIMBABWE                     9.9072e-01         8.6800e-01     4.3400e-01
#> MALAWI                       5.2558e-01         4.5600e-01     2.2800e-01
#> MOZAMBIQUE                   9.2950e-01         9.4400e-01     4.7200e-01
#> SOUTH AFRICA                 3.2564e-01         3.5200e-01     1.7600e-01
#> LESOTHO                      7.0553e-01         4.6800e-01     2.2400e-01
#> BOTSWANA                     9.8577e-01         8.4400e-01     4.2200e-01
#> SWAZILAND                    8.9999e-01         8.8400e-01     4.4200e-01
#> MOROCCO                      7.7749e-01         9.7600e-01     4.8800e-01
#> ALGERIA                      7.0057e-01         8.0400e-01     4.0200e-01
#> TUNISIA                      4.8832e-01         4.2800e-01     2.1600e-01
#> LIBYA                        2.2482e-02         4.8000e-02     2.4000e-02
#> SUDAN                        3.5067e-04         8.0000e-03     4.0000e-03
#> EGYPT                        3.1358e-04         1.6000e-02     8.0000e-03
#>                             Skewness Kurtosis
#> THE GAMBIA               -1.7643e+00   3.3422
#> MALI                     -8.2875e-01   0.8397
#> SENEGAL                  -7.1779e-01   0.3689
#> BENIN                    -8.5318e-01   0.5520
#> MAURITANIA               -9.6301e-01   1.3949
#> NIGER                    -6.3767e-01   0.2577
#> IVORY COAST              -8.3750e-01   0.6786
#> GUINEA                   -8.8298e-01   0.7691
#> BURKINA FASO             -7.7665e-01   0.5911
#> LIBERIA                  -9.7191e-01   1.2327
#> SIERRA LEONE             -1.0737e+00   0.6367
#> GHANA                    -9.6554e-01   0.8587
#> TOGO                     -1.0772e+00   1.2283
#> CAMEROON                 -8.1844e-01   1.1746
#> NIGERIA                  -6.3061e-01   0.1168
#> GABON                    -1.1947e+00   1.1977
#> CENTRAL AFRICAN REPUBLIC -6.7964e-01  -0.1428
#> CHAD                     -9.5123e-01   0.8025
#> CONGO                    -9.5359e-01   0.9288
#> ZAIRE                     6.6929e-01   0.5829
#> ANGOLA                    1.0230e+00   1.0841
#> UGANDA                    9.6607e-01   1.3331
#> KENYA                     8.5205e-01   0.5240
#> TANZANIA                  8.5767e-01   1.2127
#> BURUNDI                  -1.0634e+00   0.9779
#> RWANDA                   -9.5905e-01   0.8365
#> SOMALIA                   1.5185e+00   2.7876
#> ETHIOPIA                  9.9587e-01   0.7374
#> ZAMBIA                    4.9495e-01   0.1900
#> ZIMBABWE                 -8.0240e-01   0.4192
#> MALAWI                   -1.1246e+00   1.7228
#> MOZAMBIQUE               -7.6989e-01   0.5793
#> SOUTH AFRICA              5.6394e-01   0.1988
#> LESOTHO                  -1.7905e+00   2.9865
#> BOTSWANA                 -9.4637e-01   0.5722
#> SWAZILAND                -1.3081e+00   1.7612
#> MOROCCO                   1.2705e+00   1.5574
#> ALGERIA                   6.3958e-01   0.2677
#> TUNISIA                   1.2498e+00   1.2249
#> LIBYA                     7.3453e-01   0.5219
#> SUDAN                     4.3993e-01  -0.1019
#> EGYPT                     1.1875e+00   1.6480
```
