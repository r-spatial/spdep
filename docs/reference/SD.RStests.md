# Rao's score and adjusted Rao's score tests of linear hypotheses for spatial Durbin and spatial Durbin error models

Rao's score and adjusted Rao's score tests of linear hypotheses applied
to a fitted linear model to examine whether either the spatially lagged
dependent variable `lag` or the spatially lagged independent variable(s)
`WX` should be included in the model, or both (SDM). Adjusted tests are
provided for `lag` and `WX` adapting to the presence of the other, and a
joint test for both. The joint test is equal to the unadjusted of one
plus the adjusted of the other. In addition, draft tests are added from
Koley (2024, section 6) for spatial Durbin error models to examine
whether either the spatially lagged error `err` or the spatially lagged
independent variable(s) `WX` should be included in the model, or both
(SDEM); because of orthogonality, no adjusted tests are required.

Because the spatial lags of categorical variables are poorly understood,
from version 1.3-11 warnings are given when categorical variables are
included in the Durbin term. A discussion can be found at
<https://github.com/rsbivand/eqc25_talk>.

## Usage

``` r
SD.RStests(model, listw, zero.policy = attr(listw, "zero.policy"), test = "SDM",
 Durbin = TRUE, data = NULL)
have_factor_preds_mf(mf)
warn_factor_preds(x)
```

## Arguments

- model:

  an object of class `lm` returned by `lm`

- listw:

  a `listw` object created for example by `nb2listw`, expected to be
  row-standardised (W-style)

- zero.policy:

  default `attr(listw, "zero.policy")` as set when `listw` was created,
  if attribute not set, use global option value; if TRUE assign zero to
  the lagged value of zones without neighbours, if FALSE assign NA

- test:

  test=“SDM” computes the SDM tests, a character vector of tests
  requested chosen from SDM_RSlag, SDM_adjRSlag, SDM_RSWX, SDM_adjRSWX,
  SDM_Joint, test=“SDEM” computes the SDEM tests, a character vector of
  tests requested chosen from SDEM_RSerr, SDEM_RSWX, SDEM_Joint;
  test=“all” computes all the tests

- Durbin:

  default TRUE for Durbin models including WX; if TRUE, full spatial
  Durbin model; if a formula object, the subset of explanatory variables
  to lag

- data:

  original data object from lm() call, required if `Durbin` is a formula
  object, because the model frame must be regenerated from the original
  input data

- mf:

  `model.frame` object

- x:

  object created by `have_factor_preds_mf`

## Value

A list of class `LMtestlist` of `htest` objects, each with:

- statistic:

  the value of the Lagrange Multiplier test.

- parameter:

  number of degrees of freedom

- p.value:

  the p-value of the test.

- method:

  a character string giving the method used.

- data.name:

  a character string giving the name(s) of the data.

## References

Malabika Koley and Anil K. Bera (2024) To use, or not to use the spatial
Durbin model? – that is the question, Spatial Economic Analysis, 19:1,
30-56,
[doi:10.1080/17421772.2023.2256810](https://doi.org/10.1080/17421772.2023.2256810)
; Malabika Koley (2024) Specification testing of spatial econometric
models, Ph.D. dissertation, University of Illinois at Urbana-Champaign,
pp. 78-118, 148-156; <https://hdl.handle.net/2142/124329>.

## Author

Roger Bivand <Roger.Bivand@nhh.no>, Malabika Koley and Anil K. Bera

## Note

The results in the example below agree with those in Table 3, p. 22 in
Koley and Bera (2024).

## See also

[`lm`](https://rdrr.io/r/stats/lm.html),
[`lm.RStests`](https://r-spatial.github.io/spdep/reference/lm.RStests.md)

## Examples

``` r
columbus <- sf::st_read(system.file("shapes/columbus.gpkg", package="spData")[1])
#> Reading layer `columbus' from data source 
#>   `/home/rsb/lib/r_libs/spData/shapes/columbus.gpkg' using driver `GPKG'
#> Simple feature collection with 49 features and 20 fields
#> Geometry type: POLYGON
#> Dimension:     XY
#> Bounding box:  xmin: 5.874907 ymin: 10.78863 xmax: 11.28742 ymax: 14.74245
#> Projected CRS: Undefined Cartesian SRS with unknown unit
col.gal.nb <- read.gal(system.file("weights/columbus.gal", package="spData")[1])
col.listw <- nb2listw(col.gal.nb, style="W")
columbus$fEW <- factor(columbus$EW)
columbus$fDISCBD <- ordered(cut(columbus$DISCBD, c(0, 1.5, 3, 4.5, 6)))
f <- formula(log(CRIME) ~ INC + HOVAL + fDISCBD + fEW)
lm_obj <- lm(f, data=columbus)
summary(lm.RStests(lm_obj, col.listw, test="all"))
#>  Rao's score (a.k.a Lagrange multiplier) diagnostics for spatial
#>  dependence
#> data:  
#> model: lm(formula = f, data = columbus)
#> test weights: col.listw
#>  
#>          statistic parameter p.value
#> RSerr     1.915043         1  0.1664
#> RSlag     1.827463         1  0.1764
#> adjRSerr  0.143411         1  0.7049
#> adjRSlag  0.055831         1  0.8132
#> SARMA     1.970874         2  0.3733
res <- SD.RStests(lm_obj, col.listw, test="SDM")
#> Warning: use of spatially lagged factors (categorical variables)
#> fDISCBD, fEW
#> is not well-understood
#> Warning: In addition variable:
#> fDISCBD
#> is ordered (ordinal) with polynomial contrasts.
summary(res)
#>  Rao's score test spatial Durbin diagnostics
#> data:  
#> model: lm(formula = f, data = columbus)
#> weights: col.listw
#>  
#>              statistic parameter p.value
#> SDM_RSlag       1.8275         1  0.1764
#> SDM_adjRSlag    1.9150         1  0.1664
#> SDM_RSWX        3.6798         6  0.7199
#> SDM_adjRSWX     3.7673         6  0.7081
#> SDM_Joint       5.5948         7  0.5878
all.equal(unname(res$SDM_Joint$statistic),
 unname(res$SDM_RSlag$statistic + res$SDM_adjRSWX$statistic))
#> [1] TRUE
all.equal(unname(res$SDM_Joint$statistic),
 unname(res$SDM_adjRSlag$statistic + res$SDM_RSWX$statistic))
#> [1] TRUE
res <- SD.RStests(lm_obj, col.listw, test="SDEM")
#> Warning: use of spatially lagged factors (categorical variables)
#> fDISCBD, fEW
#> is not well-understood
#> Warning: In addition variable:
#> fDISCBD
#> is ordered (ordinal) with polynomial contrasts.
summary(res)
#>  Rao's score test spatial Durbin diagnostics
#> data:  
#> model: lm(formula = f, data = columbus)
#> weights: col.listw
#>  
#>            statistic parameter p.value
#> SDEM_RSerr    1.9150         1  0.1664
#> SDEM_RSWX     3.6798         6  0.7199
#> SDEM_Joint    5.5948         7  0.5878
all.equal(unname(res$SDEM_Joint$statistic),
 unname(res$SDEM_RSerr$statistic + res$SDEM_RSWX$statistic))
#> [1] TRUE
summary(SD.RStests(lm_obj, nb2listw(col.gal.nb, style="C"), test="all"))
#> Warning: Spatial weights matrix not row standardized
#> Warning: use of spatially lagged factors (categorical variables)
#> fDISCBD, fEW
#> is not well-understood
#> Warning: In addition variable:
#> fDISCBD
#> is ordered (ordinal) with polynomial contrasts.
#>  Rao's score test spatial Durbin diagnostics
#> data:  
#> model: lm(formula = f, data = columbus)
#> weights: nb2listw(col.gal.nb, style = "C")
#>  
#>              statistic parameter p.value  
#> SDM_RSlag       3.0731         1 0.07960 .
#> SDM_adjRSlag    3.0436         1 0.08106 .
#> SDM_RSWX        7.0092         7 0.42792  
#> SDM_adjRSWX     6.9797         7 0.43100  
#> SDM_Joint      10.0528         8 0.26134  
#> SDEM_RSerr      3.0436         1 0.08106 .
#> SDEM_RSWX       7.0092         7 0.42792  
#> SDEM_Joint     10.0528         8 0.26134  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
summary(SD.RStests(lm_obj, col.listw, test="all", Durbin= ~ INC, data=columbus))
#>  Rao's score test spatial Durbin diagnostics
#> data:  
#> model: lm(formula = f, data = columbus)
#> weights: col.listw
#> Durbin: ~ NA
#>  
#>              statistic parameter p.value
#> SDM_RSlag     1.827463         1  0.1764
#> SDM_adjRSlag  1.816700         1  0.1777
#> SDM_RSWX      0.048043         1  0.8265
#> SDM_adjRSWX   0.037280         1  0.8469
#> SDM_Joint     1.864743         2  0.3936
#> SDEM_RSerr    1.915043         1  0.1664
#> SDEM_RSWX     0.048043         1  0.8265
#> SDEM_Joint    1.963086         2  0.3747
lm_obj0 <- lm(I(scale(CRIME)) ~ 0 + I(scale(INC)) + I(scale(HOVAL)),
 data=columbus)
summary(SD.RStests(lm_obj0, col.listw, test="all"))
#>  Rao's score test spatial Durbin diagnostics
#> data:  
#> model: lm(formula = I(scale(CRIME)) ~ 0 + I(scale(INC)) +
#> I(scale(HOVAL)), data = columbus)
#> weights: col.listw
#>  
#>              statistic parameter  p.value   
#> SDM_RSlag       7.8250         1 0.005153 **
#> SDM_adjRSlag    4.6111         1 0.031765 * 
#> SDM_RSWX        6.0609         2 0.048295 * 
#> SDM_adjRSWX     2.8470         2 0.240873   
#> SDM_Joint      10.6720         3 0.013638 * 
#> SDEM_RSerr      4.6111         1 0.031765 * 
#> SDEM_RSWX       6.0609         2 0.048295 * 
#> SDEM_Joint     10.6720         3 0.013638 * 
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
columbusNA <- columbus
columbusNA$HOVAL[15] <- NA
lm_objNA <- lm(CRIME ~ INC + HOVAL, data=columbusNA)
summary(SD.RStests(lm_objNA, col.listw, test="all"))
#>  Rao's score test spatial Durbin diagnostics
#> data:  
#> model: lm(formula = CRIME ~ INC + HOVAL, data = columbusNA)
#> weights: col.listw
#>  
#>              statistic parameter  p.value   
#> SDM_RSlag       7.6010         1 0.005834 **
#> SDM_adjRSlag    4.4716         1 0.034463 * 
#> SDM_RSWX        6.0812         2 0.047807 * 
#> SDM_adjRSWX     2.9518         2 0.228576   
#> SDM_Joint      10.5527         3 0.014407 * 
#> SDEM_RSerr      4.4716         1 0.034463 * 
#> SDEM_RSWX       6.0812         2 0.047807 * 
#> SDEM_Joint     10.5527         3 0.014407 * 
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
