# Saddlepoint approximation of local Moran's Ii tests

The function implements Tiefelsdorf's application of the Saddlepoint
approximation to local Moran's Ii's reference distribution. If the model
object is of class "lm", global independence is assumed; if of class
"sarlm", global dependence is assumed to be represented by the spatial
parameter of that model. Tests are reported separately for each zone
selected, and may be summarised using `summary.localmoransad`. Values of
local Moran's Ii agree with those from
[`localmoran()`](https://r-spatial.github.io/spdep/reference/localmoran.md),
but in that function, the standard deviate - here the Saddlepoint
approximation - is based on the randomisation assumption.

## Usage

    localmoran.sad(model, select, nb, glist=NULL, style="W",
     zero.policy=NULL, alternative="two.sided", spChk=NULL,
     resfun=weighted.residuals, save.Vi=FALSE,
     tol = .Machine$double.eps^0.5, maxiter = 1000, tol.bounds=0.0001,
     save.M=FALSE, Omega = NULL)
    <!-- %as.data.frame.localmoransad(x, row.names=NULL, optional=FALSE)  -->
    # S3 method for class 'localmoransad'
    print(x, ...)
    # S3 method for class 'localmoransad'
    summary(object, ...)
    # S3 method for class 'summary.localmoransad'
    print(x, ...)
    listw2star(listw, ireg, style, n, D, a, zero.policy=attr(listw, "zero.policy"))

## Arguments

- model:

  an object of class `lm` returned by `lm` (assuming no global spatial
  autocorrelation), or an object of class `sarlm` returned by a spatial
  simultaneous autoregressive model fit (assuming global spatial
  autocorrelation represented by the model spatial coefficient); weights
  may be specified in the `lm` fit, but offsets should not be used

- select:

  an integer vector of the id. numbers of zones to be tested; if
  missing, all zones

- nb:

  a list of neighbours of class `nb`

- glist:

  a list of general weights corresponding to neighbours

- style:

  can take values W, B, C, and S

- zero.policy:

  default `attr(listw, "zero.policy")` as set when `listw` was created,
  if attribute not set, use global option value; if TRUE assign zero to
  the lagged value of zones without neighbours, if FALSE assign NA

- alternative:

  a character string specifying the alternative hypothesis, must be one
  of greater (default), less or two.sided.

- spChk:

  should the data vector names be checked against the spatial objects
  for identity integrity, TRUE, or FALSE, default NULL to use
  [`get.spChkOption()`](https://r-spatial.github.io/spdep/reference/set.spChkOption.md)

- resfun:

  default: weighted.residuals; the function to be used to extract
  residuals from the `lm` object, may be `residuals`,
  `weighted.residuals`, `rstandard`, or `rstudent`

- save.Vi:

  if TRUE, return the star-shaped weights lists for each zone tested

- tol:

  the desired accuracy (convergence tolerance) for `uniroot`

- maxiter:

  the maximum number of iterations for `uniroot`

- tol.bounds:

  offset from bounds for `uniroot`

- save.M:

  if TRUE, save a list of left and right M products in a list for the
  conditional tests, or a list of the regression model matrix components

- Omega:

  A SAR process matrix may be passed in to test an alternative
  hypothesis, for example
  `Omega <- invIrW(listw, rho=0.1); Omega <- tcrossprod(Omega)`,
  [`chol()`](https://rdrr.io/pkg/Matrix/man/chol-methods.html) is taken
  internally

- x:

  object to be printed

&nbsp;

- object:

  object to be summarised

- ...:

  arguments to be passed through

- listw:

  a `listw` object created for example by `nb2listw`

- ireg:

  a zone number

- n:

  internal value depending on listw and style

- D:

  internal value depending on listw and style

- a:

  internal value depending on listw and style

## Details

The function implements the analytical eigenvalue calculation together
with trace shortcuts given or suggested in Tiefelsdorf (2002), partly
following remarks by J. Keith Ord, and uses the Saddlepoint analytical
solution from Tiefelsdorf's SPSS code.

If a histogram of the probability values of the saddlepoint estimate for
the assumption of global independence is not approximately flat, the
assumption is probably unjustified, and re-estimation with global
dependence is recommended.

No n by n matrices are needed at any point for the test assuming no
global dependence, the star-shaped weights matrices being handled as
listw lists. When the test is made on residuals from a spatial
regression, taking a global process into account. n by n matrices are
necessary, and memory constraints may be reached for large lattices.

## Value

A list with class `localmoransad` containing "select" lists, each with
class `moransad` with the following components:

- statistic:

  the value of the saddlepoint approximation of the standard deviate of
  local Moran's Ii.

- p.value:

  the p-value of the test.

- estimate:

  the value of the observed local Moran's Ii.

- alternative:

  a character string describing the alternative hypothesis.

- method:

  a character string giving the method used.

- data.name:

  a character string giving the name(s) of the data.

- internal1:

  Saddlepoint omega, r and u

- df:

  degrees of freedom

- tau:

  maximum and minimum analytical eigenvalues

- i:

  zone tested

## References

Tiefelsdorf, M. 2002 The Saddlepoint approximation of Moran's I and
local Moran's Ii reference distributions and their numerical evaluation.
Geographical Analysis, 34, pp. 187–206.

## Author

Roger Bivand <Roger.Bivand@nhh.no>

## See also

[`localmoran`](https://r-spatial.github.io/spdep/reference/localmoran.md),
[`lm.morantest`](https://r-spatial.github.io/spdep/reference/lm.morantest.md),
[`lm.morantest.sad`](https://r-spatial.github.io/spdep/reference/lm.morantest.sad.md),
[`errorsarlm`](https://r-spatial.github.io/spdep/reference/spdep-defunct.md)

## Examples

``` r
eire <- st_read(system.file("shapes/eire.gpkg", package="spData")[1])
#> Reading layer `eire' from data source 
#>   `/home/rsb/lib/r_libs/spData/shapes/eire.gpkg' using driver `GPKG'
#> Simple feature collection with 26 features and 10 fields
#> Geometry type: MULTIPOLYGON
#> Dimension:     XY
#> Bounding box:  xmin: -4.12 ymin: 5768 xmax: 300.82 ymax: 6119.25
#> Projected CRS: Undefined Cartesian SRS with unknown unit
row.names(eire) <- as.character(eire$names)
eire.nb <- poly2nb(eire)
lw <- nb2listw(eire.nb)
e.lm <- lm(OWNCONS ~ ROADACC, data=eire)
e.locmor <- summary(localmoran.sad(e.lm, nb=eire.nb))
e.locmor
#>              Local Morans I Stand. dev. (N)      Pr. (N) Saddlepoint
#> 1 Carlow         0.21699668      0.74177148 4.582258e-01  0.95074844
#> 2 Cavan         -0.37257361     -0.81297374 1.583767e+00 -1.00603119
#> 3 Clare          0.23197510      0.49502499 6.205825e-01  0.67166518
#> 4 Cork           0.78193548      1.75999915 7.840795e-02  1.74761575
#> 5 Donegal       -1.69064059     -1.98244914 1.952571e+00 -1.72031078
#> 6 Dublin        -0.16069692     -0.15041940 1.119566e+00 -0.35212627
#> 7 Galway         1.31371473      3.34305297 8.286208e-04  2.66849536
#> 8 Kerry          0.36534866      0.58147812 5.609183e-01  0.78073279
#> 9 Kildare       -0.02557544      0.15146558 8.796085e-01  0.04167665
#> 10 Kilkenny      0.57684331      1.62868431 1.033799e-01  1.70897697
#> 11 Laoghis      -0.05951798      0.01724651 9.862400e-01 -0.12155465
#> 12 Leitrim       0.38484587      1.27548777 2.021367e-01  1.47227033
#> 13 Limerick      0.11817987      0.34038037 7.335701e-01  0.45727712
#> 14 Longford      1.41643200      3.10732224 1.887905e-03  2.51113769
#> 15 Louth         0.56242920      0.87665682 3.806731e-01  1.07441571
#> 16 Mayo          0.87572704      2.02375701 4.299516e-02  2.05251226
#> 17 Meath         0.00367856      0.16081712 8.722374e-01  0.12813539
#> 18 Monaghan      0.55098311      1.06684464 2.860420e-01  1.23999193
#> 19 Offaly        0.15155556      0.61933942 5.356928e-01  0.80786519
#> 20 Roscommon     2.04368839      6.66107106 2.718381e-11  4.53187292
#> 21 Sligo        -0.47579871     -0.73430274 1.537236e+00 -0.94578114
#> 22 Tipperary    -0.03454106      0.04843351 9.613708e-01 -0.06919691
#> 23 Waterford     0.85723423      1.98516133 4.712653e-02  1.91385108
#> 24 Westmeath     0.45138572      1.20305006 2.289569e-01  1.36017204
#> 25 Wexford       0.64371834      1.55468550 1.200210e-01  1.63188492
#> 26 Wicklow       0.02441950      0.21823347 8.272472e-01  0.21197000
#>                 Pr. (Sad) Expectation   Variance   Skewness Kurtosis    Minimum
#> 1 Carlow     3.417321e-01 -0.06471134 0.14423085 -0.7895263 7.059116  -5.461326
#> 2 Cavan      3.144006e-01 -0.04030838 0.16703857 -0.4622793 6.813004  -5.567356
#> 3 Clare      5.017969e-01 -0.04219377 0.30674823 -0.3579185 6.761932  -7.406885
#> 4 Cork       8.053059e-02 -0.04363826 0.22003250 -0.4363325 6.799081  -6.360927
#> 5 Donegal    8.537596e-02 -0.11935726 0.62821008 -0.7004303 6.979051 -11.236382
#> 6 Dublin     7.247436e-01 -0.08240217 0.27093033 -0.7352902 7.009201  -7.420687
#> 7 Galway     7.619183e-03 -0.04389010 0.16491503 -0.5059997 6.838308  -5.573707
#> 8 Kerry      4.349597e-01 -0.03212369 0.46724761 -0.2212594 6.714810  -8.915104
#> 9 Kildare    9.667565e-01 -0.07623733 0.11187546 -1.0416026 7.339962  -4.999621
#> 10 Kilkenny  8.745521e-02 -0.05911601 0.15247016 -0.7040737 6.982131  -5.538889
#> 11 Laoghis   9.032517e-01 -0.06620117 0.15016408 -0.7915137 7.061015  -5.574277
#> 12 Leitrim   1.409479e-01 -0.06828680 0.12621127 -0.8864130 7.157471  -5.186756
#> 13 Limerick  6.474719e-01 -0.04172917 0.22070746 -0.4167943 6.789132  -6.348870
#> 14 Longford  1.203427e-02 -0.04055162 0.21985521 -0.4059156 6.783792  -6.324457
#> 15 Louth     2.826364e-01 -0.04386717 0.47831137 -0.2983105 6.738631  -9.149779
#> 16 Mayo      4.011990e-02 -0.10340651 0.23408152 -0.9804377 7.264328  -7.165847
#> 17 Meath     8.980418e-01 -0.04799533 0.10324707 -0.6949022 6.974408  -4.551176
#> 18 Monaghan  2.149784e-01 -0.04136950 0.30828917 -0.3501006 6.758633  -7.415047
#> 19 Offaly    4.191682e-01 -0.04662212 0.10238870 -0.6782821 6.960676  -4.519986
#> 20 Roscommon 5.846302e-06 -0.04435301 0.09826301 -0.6591526 6.945294  -4.414164
#> 21 Sligo     3.442602e-01 -0.11041170 0.24760301 -1.0156782 7.307311  -7.409091
#> 22 Tipperary 9.448329e-01 -0.04872132 0.08571886 -0.7716947 7.042300  -4.198353
#> 23 Waterford 5.563919e-02 -0.05435079 0.21086415 -0.5533724 6.868342  -6.353514
#> 24 Westmeath 1.737755e-01 -0.04181671 0.16806722 -0.4779156 6.821788  -5.599610
#> 25 Wexford   1.027037e-01 -0.05432673 0.20159596 -0.5654870 6.876461  -6.225016
#> 26 Wicklow   8.321304e-01 -0.07022471 0.18808122 -0.7515626 7.023792  -6.198974
#>               Maximum        omega       sad.r       sad.u
#> 1 Carlow     3.908254  0.071122863  0.71853146  0.84900467
#> 2 Cavan      4.599954 -0.051532330 -0.69003582 -0.85816092
#> 3 Clare      6.394234  0.031939551  0.46982957  0.51656357
#> 4 Cork       5.313609  0.083281515  1.43928884  2.24323856
#> 5 Donegal    8.371808 -0.039371561 -1.39288713 -2.19776923
#> 6 Dublin     5.443035 -0.010354512 -0.14127218 -0.14554367
#> 7 Galway     4.520344  0.134873939  2.39189491  4.63522665
#> 8 Kerry      8.144135  0.028431952  0.54051534  0.61545648
#> 9 Kildare    3.169925  0.018723213  0.14926096  0.14688325
#> 10 Kilkenny  4.120105  0.105365819  1.40668493  2.15214600
#> 11 Laoghis   3.985449  0.001724096  0.01660985  0.01657177
#> 12 Leitrim   3.547872  0.109608339  1.18229964  1.66578026
#> 13 Limerick  5.347370  0.027585567  0.32820959  0.34241163
#> 14 Longford  5.351218  0.108368224  2.22687206  4.19385821
#> 15 Louth     8.096967  0.037745731  0.78968254  0.98878788
#> 16 Mayo      4.684091  0.105317102  1.75922892  2.94711915
#> 17 Meath     3.399289  0.020074178  0.15693992  0.15623206
#> 18 Monaghan  6.422179  0.053418031  0.94193428  1.24723730
#> 19 Offaly    3.401055  0.071432961  0.59871278  0.67858124
#> 20 Roscommon 3.349691  0.357114199  4.33885200 10.02517077
#> 21 Sligo     4.759210 -0.036185300 -0.60775394 -0.74635955
#> 22 Tipperary 3.029041  0.006482570  0.04683010  0.04657633
#> 23 Waterford 5.049095  0.093409563  1.61022194  2.62552879
#> 24 Westmeath 4.596009  0.080097979  1.05978192  1.45704611
#> 25 Wexford   4.921174  0.085478829  1.32647942  1.98902088
#> 26 Wicklow   4.513581  0.020437514  0.21417067  0.21406975
mean(e.locmor[,1])
#> [1] 0.3366057
sum(e.locmor[,1])/Szero(lw)
#> [1] 0.3366057
lm.morantest(e.lm, lw)
#> 
#>  Global Moran I for regression residuals
#> 
#> data:  
#> model: lm(formula = OWNCONS ~ ROADACC, data = eire)
#> weights: lw
#> 
#> Moran I statistic standard deviate = 3.2575, p-value = 0.0005619
#> alternative hypothesis: greater
#> sample estimates:
#> Observed Moran I      Expectation         Variance 
#>       0.33660565      -0.05877741       0.01473183 
#> 
# note equality for mean() only when the sum of weights equals
# the number of observations (thanks to Juergen Symanzik)
hist(e.locmor[,"Pr. (Sad)"])

e.wlm <- lm(OWNCONS ~ ROADACC, data=eire, weights=RETSALE)
e.locmorw1 <- summary(localmoran.sad(e.wlm, nb=eire.nb, resfun=weighted.residuals))
e.locmorw1
#>              Local Morans I Stand. dev. (N)      Pr. (N) Saddlepoint
#> 1 Carlow        0.160490657      0.41009538 6.817360e-01  0.57729730
#> 2 Cavan        -0.144663771     -0.27825575 1.219184e+00 -0.45305168
#> 3 Clare         0.301676582      0.59021701 5.550452e-01  0.78987327
#> 4 Cork          1.520488153      3.49016956 4.827142e-04  2.83665634
#> 5 Donegal      -0.435899664     -0.44665946 1.344879e+00 -0.67136090
#> 6 Dublin       -0.340222996     -0.60886675 1.457387e+00 -0.85658016
#> 7 Galway        1.469009037      3.64546192 2.669119e-04  2.82274510
#> 8 Kerry         0.268957279      0.49517211 6.204786e-01  0.66861012
#> 9 Kildare       0.002996577      0.12659200 8.992633e-01  0.09978249
#> 10 Kilkenny     0.497968801      1.20610829 2.277757e-01  1.34091379
#> 11 Laoghis     -0.026328708     -0.01713904 1.013674e+00 -0.06511750
#> 12 Leitrim      0.277306386      0.80644537 4.199861e-01  1.01176491
#> 13 Limerick     0.263122939      0.67591893 4.990921e-01  0.87759628
#> 14 Longford     0.605211711      1.23874274 2.154408e-01  1.35375850
#> 15 Louth        0.359469683      0.53572990 5.921452e-01  0.73280188
#> 16 Mayo         1.742480472      3.68118786 2.321499e-04  3.01500163
#> 17 Meath       -0.112171618     -0.26680698 1.210382e+00 -0.45470531
#> 18 Monaghan     0.241044892      0.46013751 6.454175e-01  0.64094159
#> 19 Offaly       0.129191178      0.40988983 6.818868e-01  0.57205532
#> 20 Roscommon    1.469671750      4.39927211 1.086146e-05  3.13215267
#> 21 Sligo        0.095559559      0.31036587 7.562827e-01  0.38503064
#> 22 Tipperary    0.088352646      0.39958273 6.894639e-01  0.52592105
#> 23 Waterford    1.212658861      2.69846035 6.966104e-03  2.29292584
#> 24 Westmeath    0.256165310      0.61222265 5.403905e-01  0.81620859
#> 25 Wexford      0.619197906      1.31967483 1.869436e-01  1.42553832
#> 26 Wicklow      0.006531073      0.13142161 8.954418e-01  0.10986725
#>                Pr. (Sad) Expectation   Variance   Skewness Kurtosis    Minimum
#> 1 Carlow     0.563738639 -0.01840558 0.19029731 -0.1986977 6.709176  -5.665283
#> 2 Cavan      0.650511526 -0.02184465 0.19482455 -0.2329761 6.717975  -5.769371
#> 3 Clare      0.429601786 -0.02782578 0.31166901 -0.2346279 6.718435  -7.299407
#> 4 Cork       0.004558865 -0.07147410 0.20805228 -0.7280203 7.002789  -6.495227
#> 5 Donegal    0.501990650 -0.07638800 0.64784767 -0.4450323 6.803659 -10.931345
#> 6 Dublin     0.391676955 -0.11187469 0.14065367 -1.3324828 7.767845  -5.846924
#> 7 Galway     0.004761441 -0.04407925 0.17227545 -0.4973247 6.833103  -5.688263
#> 8 Kerry      0.503744214 -0.05644300 0.43184063 -0.4031537 6.782459  -8.859407
#> 9 Kildare    0.920517009 -0.04145378 0.12329287 -0.5519841 6.867423  -4.857131
#> 10 Kilkenny  0.179948446 -0.02321457 0.18672759 -0.2528322 6.723714  -5.669156
#> 11 Laoghis   0.948080445 -0.01882531 0.19166504 -0.2024943 6.710081  -5.689691
#> 12 Leitrim   0.311650473 -0.03314655 0.14819776 -0.4041387 6.782933  -5.190860
#> 13 Limerick  0.380162847 -0.04800976 0.21188592 -0.4885390 6.827925  -6.298874
#> 14 Longford  0.175813436 -0.01321709 0.24923967 -0.1247567 6.694960  -6.392250
#> 15 Louth     0.463679274 -0.02272663 0.50895713 -0.1500895 6.699098  -9.179410
#> 16 Mayo      0.002569779 -0.09458039 0.24904116 -0.8745463 7.144787  -7.272523
#> 17 Meath     0.649321260 -0.02548973 0.10555100 -0.3685280 6.766525  -4.353122
#> 18 Monaghan  0.521560632 -0.02333303 0.33012248 -0.1912612 6.707452  -7.451275
#> 19 Offaly    0.567284497 -0.01927989 0.13120470 -0.2505068 6.723017  -4.750083
#> 20 Roscommon 0.001735296 -0.02385522 0.11525625 -0.3302939 6.750604  -4.517676
#> 21 Sligo     0.700214717 -0.06927688 0.28207107 -0.6087587 6.906925  -7.416881
#> 22 Tipperary 0.598943063 -0.03719246 0.09871568 -0.5534442 6.868390  -4.347214
#> 23 Waterford 0.021852275 -0.04027429 0.21558757 -0.4070996 6.784366  -6.264083
#> 24 Westmeath 0.414380842 -0.01522814 0.19650786 -0.1618342 6.701275  -5.716689
#> 25 Wexford   0.154001651 -0.02786850 0.24041626 -0.2674357 6.728234  -6.450174
#> 26 Wicklow   0.912514660 -0.04751424 0.16911517 -0.5404006 6.859847  -5.677417
#>               Maximum        omega       sad.r       sad.u
#> 1 Carlow     5.223549  0.033637462  0.38794383  0.41751442
#> 2 Cavan      5.245099 -0.022294725 -0.26091513 -0.27432856
#> 3 Clare      6.631589  0.035302524  0.54866959  0.62630505
#> 4 Cork       4.779848  0.134384022  2.56890121  5.11058941
#> 5 Donegal    9.098033 -0.017507303 -0.40402088 -0.45010348
#> 6 Dublin     3.161931 -0.040828131 -0.50710580 -0.60543071
#> 7 Galway     4.630361  0.138454493  2.55314856  5.08170760
#> 8 Kerry      7.504775  0.027195108  0.47167584  0.51758880
#> 9 Kildare    3.862240  0.014230854  0.12280697  0.12246022
#> 10 Kilkenny  5.112007  0.071351828  1.03314661  1.41988761
#> 11 Laoghis   5.237883 -0.001502568 -0.01645549 -0.01646868
#> 12 Leitrim   4.395343  0.066207108  0.74158109  0.90609823
#> 13 Limerick  5.146640  0.050265468  0.63703133  0.74253221
#> 14 Longford  6.075040  0.060534639  1.04112815  1.44165700
#> 15 Louth     8.633971  0.025231866  0.49793592  0.55971045
#> 16 Mayo      5.002594  0.135902733  2.75603129  5.62669411
#> 17 Meath     3.741368 -0.028724497 -0.24912602 -0.26221742
#> 18 Monaghan  6.891283  0.028042744  0.43280552  0.47360371
#> 19 Offaly    4.287365  0.040883462  0.38909358  0.41780262
#> 20 Roscommon 3.945151  0.180270100  2.87643937  6.00202551
#> 21 Sligo     5.754236  0.023145822  0.30299868  0.31062426
#> 22 Tipperary 3.454595  0.048932861  0.38768333  0.40902706
#> 23 Waterford 5.297500  0.101874758  1.99944904  3.59542986
#> 24 Westmeath 5.351214  0.044919886  0.56403221  0.65024353
#> 25 Wexford   5.781330  0.066083430  1.11524637  1.57637406
#> 26 Wicklow   4.537075  0.012609906  0.12748740  0.12720134
e.locmorw2 <- summary(localmoran.sad(e.wlm, nb=eire.nb, resfun=rstudent))
e.locmorw2
#>              Local Morans I Stand. dev. (N)      Pr. (N) Saddlepoint
#> 1 Carlow        0.121997164      0.32185426 0.7475631028  0.45871930
#> 2 Cavan        -0.110588620     -0.20105599 1.1593452023 -0.34808741
#> 3 Clare         0.308759283      0.60290381 0.5465726820  0.80379674
#> 4 Cork          1.404203115      3.23522979 0.0012154487  2.69598105
#> 5 Donegal      -0.421182109     -0.42837428 1.3316213558 -0.65187732
#> 6 Dublin       -0.443534307     -0.88433550 1.6234848706 -1.05167518
#> 7 Galway        1.376190254      3.42183495 0.0006220006  2.70663064
#> 8 Kerry         0.241623387      0.45357724 0.6501331373  0.61515126
#> 9 Kildare       0.008842913      0.14324202 0.8860990418  0.12690939
#> 10 Kilkenny     0.382893124      0.93980330 0.3473184643  1.12810394
#> 11 Laoghis     -0.019904450     -0.00246494 1.0019667357 -0.04176107
#> 12 Leitrim      0.200483743      0.60688769 0.5439254449  0.80241456
#> 13 Limerick     0.238810092      0.62310059 0.5332184222  0.81824678
#> 14 Longford     0.466942853      0.96178354 0.3361583457  1.14088225
#> 15 Louth        0.271710988      0.41271732 0.6798137369  0.58520655
#> 16 Mayo         1.705850161      3.60778634 0.0003088206  2.97346461
#> 17 Meath       -0.160603054     -0.41587901 1.3225014587 -0.63254531
#> 18 Monaghan     0.181378604      0.35629110 0.7216225742  0.50697175
#> 19 Offaly       0.107919508      0.35116431 0.7254650854  0.49312417
#> 20 Roscommon    1.249826037      3.75170355 0.0001756371  2.81402777
#> 21 Sligo        0.079431883      0.27999953 0.7794778655  0.33741411
#> 22 Tipperary    0.085153638      0.38940098 0.6969795442  0.51126046
#> 23 Waterford    1.022922345      2.28982205 0.0220316342  2.06385578
#> 24 Westmeath    0.201539331      0.48899468 0.6248454598  0.67791734
#> 25 Wexford      0.481114325      1.03805702 0.2992435100  1.21094348
#> 26 Wicklow     -0.011694661      0.08710223 0.9305902563  0.03786608
#>                Pr. (Sad) Expectation   Variance   Skewness Kurtosis    Minimum
#> 1 Carlow     0.646435751 -0.01840558 0.19029731 -0.1986977 6.709176  -5.665283
#> 2 Cavan      0.727774542 -0.02184465 0.19482455 -0.2329761 6.717975  -5.769371
#> 3 Clare      0.421514374 -0.02782578 0.31166901 -0.2346279 6.718435  -7.299407
#> 4 Cork       0.007018166 -0.07147410 0.20805228 -0.7280203 7.002789  -6.495227
#> 5 Donegal    0.514480318 -0.07638800 0.64784767 -0.4450323 6.803659 -10.931345
#> 6 Dublin     0.292948602 -0.11187469 0.14065367 -1.3324828 7.767845  -5.846924
#> 7 Galway     0.006796983 -0.04407925 0.17227545 -0.4973247 6.833103  -5.688263
#> 8 Kerry      0.538454834 -0.05644300 0.43184063 -0.4031537 6.782459  -8.859407
#> 9 Kildare    0.899012117 -0.04145378 0.12329287 -0.5519841 6.867423  -4.857131
#> 10 Kilkenny  0.259276027 -0.02321457 0.18672759 -0.2528322 6.723714  -5.669156
#> 11 Laoghis   0.966689170 -0.01882531 0.19166504 -0.2024943 6.710081  -5.689691
#> 12 Leitrim   0.422313193 -0.03314655 0.14819776 -0.4041387 6.782933  -5.190860
#> 13 Limerick  0.413216288 -0.04800976 0.21188592 -0.4885390 6.827925  -6.298874
#> 14 Longford  0.253918927 -0.01321709 0.24923967 -0.1247567 6.694960  -6.392250
#> 15 Louth     0.558408841 -0.02272663 0.50895713 -0.1500895 6.699098  -9.179410
#> 16 Mayo      0.002944584 -0.09458039 0.24904116 -0.8745463 7.144787  -7.272523
#> 17 Meath     0.527030613 -0.02548973 0.10555100 -0.3685280 6.766525  -4.353122
#> 18 Monaghan  0.612174643 -0.02333303 0.33012248 -0.1912612 6.707452  -7.451275
#> 19 Offaly    0.621924855 -0.01927989 0.13120470 -0.2505068 6.723017  -4.750083
#> 20 Roscommon 0.004892500 -0.02385522 0.11525625 -0.3302939 6.750604  -4.517676
#> 21 Sligo     0.735804747 -0.06927688 0.28207107 -0.6087587 6.906925  -7.416881
#> 22 Tipperary 0.609168686 -0.03719246 0.09871568 -0.5534442 6.868390  -4.347214
#> 23 Waterford 0.039031392 -0.04027429 0.21558757 -0.4070996 6.784366  -6.264083
#> 24 Westmeath 0.497824103 -0.01522814 0.19650786 -0.1618342 6.701275  -5.716689
#> 25 Wexford   0.225917068 -0.02786850 0.24041626 -0.2674357 6.728234  -6.450174
#> 26 Wicklow   0.969794461 -0.04751424 0.16911517 -0.5404006 6.859847  -5.677417
#>               Maximum         omega        sad.r        sad.u
#> 1 Carlow     5.223549  0.0272472954  0.306852385  0.321490350
#> 2 Cavan      5.245099 -0.0166242430 -0.190223798 -0.196022741
#> 3 Clare      6.631589  0.0358559846  0.559584155  0.641527095
#> 4 Cork       4.779848  0.1283688794  2.421993662  4.702934834
#> 5 Donegal    9.098033 -0.0169574018 -0.388662370 -0.430528174
#> 6 Dublin     3.161931 -0.0507848831 -0.696824139 -0.892299501
#> 7 Galway     4.630361  0.1333224875  2.431707338  4.745166958
#> 8 Kerry      7.504775  0.0252982005  0.433626352  0.469138111
#> 9 Kildare    3.862240  0.0161278548  0.139072871  0.138837812
#> 10 Kilkenny  5.112007  0.0622219791  0.835166590  1.066654825
#> 11 Laoghis   5.237883 -0.0002164925 -0.002368023 -0.002368244
#> 12 Leitrim   4.395343  0.0543844388  0.571625643  0.652237768
#> 13 Limerick  5.146640  0.0473787285  0.590605630  0.675595582
#> 14 Longford  6.075040  0.0528365356  0.840671447  1.082014178
#> 15 Louth     8.633971  0.0204960239  0.389078238  0.419930705
#> 16 Mayo      5.002594  0.1340509166  2.712635049  5.503894365
#> 17 Meath     3.741368 -0.0415411327 -0.379736494 -0.417998694
#> 18 Monaghan  6.891283  0.0226076552  0.338579266  0.358443998
#> 19 Offaly    4.287365  0.0357600182  0.335019324  0.353243066
#> 20 Roscommon 3.945151  0.1623235556  2.543240406  5.063795948
#> 21 Sligo     5.754236  0.0209568244  0.273376649  0.278204612
#> 22 Tipperary 3.454595  0.0478130890  0.377982061  0.397511387
#> 23 Waterford 5.297500  0.0941811771  1.762051071  2.998984773
#> 24 Westmeath 5.351214  0.0379323563  0.457473561  0.506014991
#> 25 Wexford   5.781330  0.0582538352  0.911459483  1.197527949
#> 26 Wicklow   4.537075  0.0083101099  0.084281843  0.083952776
run <- FALSE
if (requireNamespace("spatialreg", quietly=TRUE)) run <- TRUE
if (run) {
e.errorsar <- spatialreg::errorsarlm(OWNCONS ~ ROADACC, data=eire,
  listw=lw)
if (packageVersion("spatialreg") < "1.1.7")
  spatialreg::print.sarlm(e.errorsar)
else
  print(e.errorsar)
}
#> 
#> Call:
#> spatialreg::errorsarlm(formula = OWNCONS ~ ROADACC, data = eire, 
#>     listw = lw)
#> Type: error 
#> 
#> Coefficients:
#>      lambda (Intercept)     ROADACC 
#> 0.783970988 2.892719568 0.002800913 
#> 
#> Log likelihood: -64.12465 
if (run) {
lm.target <- lm(e.errorsar$tary ~ e.errorsar$tarX - 1)
Omega <- tcrossprod(spatialreg::invIrW(lw, rho=e.errorsar$lambda))
e.clocmor <- summary(localmoran.sad(lm.target, nb=eire.nb, Omega=Omega))
e.clocmor
}
#>              Local Morans I Stand. dev. (N)   Pr. (N) Saddlepoint    Pr. (Sad)
#> 1 Carlow         0.17958462      -0.5947094 1.4479623 -0.01890335 0.9849182099
#> 2 Cavan         -0.24752628      -0.9775195 1.6716880 -1.30903979 0.1905208729
#> 3 Clare         -0.27901334      -0.9987481 1.6820833 -1.41535136 0.1569655036
#> 4 Cork           0.37808655      -0.5304596 1.4042067  0.06399747 0.9489722429
#> 5 Donegal       -1.01894688      -1.4634181 1.8566470 -1.88871177 0.0589304608
#> 6 Dublin        -0.18171297      -0.9790853 1.6724621 -1.11436912 0.2651208632
#> 7 Galway         1.02193390       0.9527103 0.3407369  1.44844014 0.1474939944
#> 8 Kerry         -0.94967914      -1.1410950 1.7461696 -2.37909992 0.0173549707
#> 9 Kildare        0.07005053      -0.6531380 1.4863327 -0.10338270 0.9176592422
#> 10 Kilkenny      0.43022231      -0.3097714 1.2432652  0.46873471 0.6392592764
#> 11 Laoghis      -0.12239133      -0.9818401 1.6738214 -0.97718184 0.3284791347
#> 12 Leitrim      -0.24203970      -0.9591963 1.6625401 -1.53246533 0.1254076471
#> 13 Limerick     -0.03214546      -0.7729433 1.5604441 -0.80301568 0.4219656706
#> 14 Longford      0.38307454      -0.4899146 1.3758057  0.16621954 0.8679841830
#> 15 Louth         0.21301552      -0.6509309 1.4849089 -0.21224290 0.8319175423
#> 16 Mayo          0.93971200       0.3895188 0.6968924  1.15011408 0.2500968881
#> 17 Meath         0.12484415      -0.6182203 1.4635699 -0.09720926 0.9225602153
#> 18 Monaghan     -0.16109919      -0.8608004 1.6106520 -1.20575831 0.2279106935
#> 19 Offaly       -0.00632492      -0.7740865 1.5611204 -0.33915072 0.7344961920
#> 20 Roscommon     1.02089429       0.6454801 0.5186161  1.21215474 0.2254531514
#> 21 Sligo        -2.01629233      -1.8957526 1.9420072 -3.59571732 0.0003234989
#> 22 Tipperary    -0.10810709      -0.9082666 1.6362626 -1.09721448 0.2725476428
#> 23 Waterford     0.44099279      -0.3488980 1.2728341  0.41219202 0.6801986840
#> 24 Westmeath    -0.06329661      -0.8267635 1.5916289 -0.81444158 0.4153920275
#> 25 Wexford       0.30764883      -0.5219835 1.3983182  0.12395154 0.9013536459
#> 26 Wicklow      -0.01696406      -0.7708237 1.5591886 -0.68244661 0.4949565909
#>              Expectation   Variance Skewness  Kurtosis       Minimum   Maximum
#> 1 Carlow       0.8594106  1.3067303 2.362830 10.416149  8.621901e-17 20.625856
#> 2 Cavan        1.1247112  1.9706416 2.347202 10.347468  1.645654e+00 25.347414
#> 3 Clare        1.2679076  2.3989670 2.334022 10.287656  2.485514e+00 27.944268
#> 4 Cork         1.2841620  2.9175890 2.362830 10.416149 -1.812003e-16 30.819888
#> 5 Donegal      1.7253177  3.5165365 2.172101  9.492322  8.227735e+00 33.179890
#> 6 Dublin       0.9462667  1.3272768 2.331563 10.276363  1.928660e+00 20.781740
#> 7 Galway       0.4507425  0.3594525 2.362830 10.416149  5.907989e-16 10.817821
#> 8 Kerry        2.6933183 10.1923447 2.307813 10.165571  7.170614e+00 57.469026
#> 9 Kildare      0.5337376  0.5040110 2.362830 10.416149 -1.337526e-16 12.809701
#> 10 Kilkenny    0.7317135  0.9472545 2.362830 10.416149  1.706968e-15 17.561125
#> 11 Laoghis     0.6035313  0.5466372 2.336004 10.296727  1.143655e+00 13.341096
#> 12 Leitrim     1.2258674  2.3419747 2.347302 10.347916  1.788145e+00 27.632674
#> 13 Limerick    1.5147254  4.0050955 2.362668 10.415466  2.341575e-01 36.119253
#> 14 Longford    1.0996747  2.1395031 2.362830 10.416149  1.538224e-16 26.392192
#> 15 Louth       1.5875248  4.4588771 2.362830 10.416149 -9.743415e-17 38.100594
#> 16 Mayo        0.6190019  0.6779043 2.362830 10.416149 -8.612978e-16 14.856045
#> 17 Meath       0.7025948  0.8733620 2.362830 10.416149  1.089845e-16 16.862274
#> 18 Monaghan    1.5519577  3.9603958 2.357838 10.394611  1.306155e+00 35.940829
#> 19 Offaly      0.3048481  0.1615939 2.362560 10.415009  6.073803e-02  7.255616
#> 20 Roscommon   0.5492905  0.5338125 2.362830 10.416149  7.680029e-16 13.182973
#> 21 Sligo       2.3130311  5.2152695 1.886102  7.938452  1.739456e+01 38.118185
#> 22 Tipperary   0.7147979  0.8208668 2.354040 10.377921  7.922254e-01 16.362924
#> 23 Waterford   0.8228671  1.1979644 2.362830 10.416149  1.480772e-16 19.748811
#> 24 Westmeath   0.8426747  1.2007879 2.360962 10.408167  4.379239e-01 19.786270
#> 25 Wexford     1.0063848  1.7918953 2.362830 10.416149  6.041304e-16 24.153236
#> 26 Wicklow     0.9927147  1.7157589 2.362598 10.415171  1.833509e-01 23.641801
#>                     omega       sad.r       sad.u
#> 1 Carlow     -0.017838447 -0.38297869 -0.33313448
#> 2 Cavan      -0.066132684 -1.21593175 -1.36168605
#> 3 Clare      -0.066612498 -1.30513958 -1.50704646
#> 4 Cork       -0.010368850 -0.33915229 -0.29581081
#> 5 Donegal    -0.045800559 -1.66568212 -2.41507712
#> 6 Dublin     -0.056683668 -1.03296356 -1.12358099
#> 7 Galway      0.031762997  1.14153962  1.62047125
#> 8 Kerry      -0.079689480 -2.19374456 -3.29442701
#> 9 Kildare    -0.025051491 -0.36678693 -0.33300845
#> 10 Kilkenny   0.003214928  0.09217308  0.09542846
#> 11 Laoghis   -0.065438986 -0.89335507 -0.96282454
#> 12 Leitrim   -0.092395573 -1.43474097 -1.65069014
#> 13 Limerick  -0.047140362 -0.97737323 -0.82423542
#> 14 Longford  -0.007176482 -0.22762630 -0.20810767
#> 15 Louth     -0.018552546 -0.58005974 -0.46861210
#> 16 Mayo       0.019992304  0.84088281  1.09059594
#> 17 Meath     -0.026958325 -0.44558661 -0.38151850
#> 18 Monaghan  -0.062388996 -1.21509194 -1.20138917
#> 19 Offaly    -0.044380046 -0.41568118 -0.40266554
#> 20 Roscommon  0.023338184  0.88169755  1.17993235
#> 21 Sligo     -0.131646728 -3.40243100 -6.56750200
#> 22 Tipperary -0.089378129 -1.06748205 -1.10190609
#> 23 Waterford  0.001151041  0.03534520  0.03581914
#> 24 Westmeath -0.058470593 -0.88397327 -0.83127666
#> 25 Wexford   -0.009550496 -0.26689929 -0.24045996
#> 26 Wicklow   -0.050401553 -0.84395069 -0.73641410
if (run) {
hist(e.clocmor[,"Pr. (Sad)"])
}
```
