# Moran's I test for residual spatial autocorrelation

Moran's I test for spatial autocorrelation in residuals from an
estimated linear model ([`lm()`](https://rdrr.io/r/stats/lm.html)).

## Usage

``` r
lm.morantest(model, listw, zero.policy=attr(listw, "zero.policy"),
 alternative = "greater", spChk=NULL, resfun=weighted.residuals, naSubset=TRUE)
```

## Arguments

- model:

  an object of class `lm` returned by `lm`; weights may be specified in
  the `lm` fit, but offsets should not be used

- listw:

  a `listw` object created for example by `nb2listw`

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

- resfun:

  default: weighted.residuals; the function to be used to extract
  residuals from the `lm` object, may be `residuals`,
  `weighted.residuals`, `rstandard`, or `rstudent`

- naSubset:

  default TRUE to subset listw object for omitted observations in model
  object (this is a change from earlier behaviour, when the
  `model$na.action` component was ignored, and the listw object had to
  be subsetted by hand)

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

  a character string giving the method used.

- data.name:

  a character string giving the name(s) of the data.

## References

Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 203,

## Author

Roger Bivand <Roger.Bivand@nhh.no>

## Details

If factors or ordered factors are among explanatory variables in the
regression, the coding of the contrasts affects the outcomes of
significance tests.

## See also

[`lm.LMtests`](https://r-spatial.github.io/spdep/reference/lm.RStests.md),
[`lm`](https://rdrr.io/r/stats/lm.html)

## Examples

``` r
data(oldcol)
oldcrime1.lm <- lm(CRIME ~ 1, data = COL.OLD)
oldcrime.lm <- lm(CRIME ~ HOVAL + INC, data = COL.OLD)
lm.morantest(oldcrime.lm, nb2listw(COL.nb, style="W"))
#> 
#>  Global Moran I for regression residuals
#> 
#> data:  
#> model: lm(formula = CRIME ~ HOVAL + INC, data = COL.OLD)
#> weights: nb2listw(COL.nb, style = "W")
#> 
#> Moran I statistic standard deviate = 2.9539, p-value = 0.001569
#> alternative hypothesis: greater
#> sample estimates:
#> Observed Moran I      Expectation         Variance 
#>      0.235638354     -0.033302866      0.008289408 
#> 
lm.LMtests(oldcrime.lm, nb2listw(COL.nb, style="W"))
#> Please update scripts to use lm.RStests in place of lm.LMtests
#> 
#>  Rao's score (a.k.a Lagrange multiplier) diagnostics for spatial
#>  dependence
#> 
#> data:  
#> model: lm(formula = CRIME ~ HOVAL + INC, data = COL.OLD)
#> test weights: listw
#> 
#> RSerr = 5.7231, df = 1, p-value = 0.01674
#> 
lm.morantest(oldcrime.lm, nb2listw(COL.nb, style="S"))
#> 
#>  Global Moran I for regression residuals
#> 
#> data:  
#> model: lm(formula = CRIME ~ HOVAL + INC, data = COL.OLD)
#> weights: nb2listw(COL.nb, style = "S")
#> 
#> Moran I statistic standard deviate = 3.1745, p-value = 0.0007504
#> alternative hypothesis: greater
#> sample estimates:
#> Observed Moran I      Expectation         Variance 
#>      0.239317561     -0.033431740      0.007381982 
#> 
lm.morantest(oldcrime1.lm, nb2listw(COL.nb, style="W"))
#> 
#>  Global Moran I for regression residuals
#> 
#> data:  
#> model: lm(formula = CRIME ~ 1, data = COL.OLD)
#> weights: nb2listw(COL.nb, style = "W")
#> 
#> Moran I statistic standard deviate = 5.6754, p-value = 6.92e-09
#> alternative hypothesis: greater
#> sample estimates:
#> Observed Moran I      Expectation         Variance 
#>      0.510951264     -0.020833333      0.008779831 
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
oldcrime.wlm <- lm(CRIME ~ HOVAL + INC, data = COL.OLD,
 weights = I(1/AREA_PL))
lm.morantest(oldcrime.wlm, nb2listw(COL.nb, style="W"),
 resfun=weighted.residuals)
#> 
#>  Global Moran I for regression residuals
#> 
#> data:  
#> model: lm(formula = CRIME ~ HOVAL + INC, data = COL.OLD, weights =
#> I(1/AREA_PL))
#> weights: nb2listw(COL.nb, style = "W")
#> 
#> Moran I statistic standard deviate = 3.0141, p-value = 0.001289
#> alternative hypothesis: greater
#> sample estimates:
#> Observed Moran I      Expectation         Variance 
#>      0.241298974     -0.032224366      0.008235091 
#> 
lm.morantest(oldcrime.wlm, nb2listw(COL.nb, style="W"),
 resfun=rstudent)
#> 
#>  Global Moran I for regression residuals
#> 
#> data:  
#> model: lm(formula = CRIME ~ HOVAL + INC, data = COL.OLD, weights =
#> I(1/AREA_PL))
#> weights: nb2listw(COL.nb, style = "W")
#> 
#> Moran I statistic standard deviate = 2.822, p-value = 0.002387
#> alternative hypothesis: greater
#> sample estimates:
#> Observed Moran I      Expectation         Variance 
#>      0.223860298     -0.032224366      0.008235091 
#> 
run <- require("codingMatrices", quietly=TRUE)
if (run) {
COL.OLD$fEW <- factor(COL.OLD$EW)
COL.OLD$fDISCBD <- ordered(cut(COL.OLD$DISCBD, c(0, 1.5, 3, 4.5, 6)))
f <- formula(CRIME ~ INC + HOVAL + fDISCBD*fEW)
lw <- nb2listw(COL.nb, style="W")
# default codings
lm.morantest(lm(f, data=COL.OLD, contrasts=list(fDISCBD="contr.poly", fEW="contr.treatment")), lw)
}
#> 
#>  Global Moran I for regression residuals
#> 
#> data:  
#> model: lm(formula = f, data = COL.OLD, contrasts = list(fDISCBD =
#> "contr.poly", fEW = "contr.treatment"))
#> weights: lw
#> 
#> Moran I statistic standard deviate = 0.98381, p-value = 0.1626
#> alternative hypothesis: greater
#> sample estimates:
#> Observed Moran I      Expectation         Variance 
#>     -0.020008220     -0.099678739      0.006558024 
#> 
if (run) {
# use codingMatrices::code_diff for ordered factor
lm.morantest(lm(f, data=COL.OLD, contrasts=list(fDISCBD="code_diff", fEW="contr.treatment")), lw)
}
#> 
#>  Global Moran I for regression residuals
#> 
#> data:  
#> model: lm(formula = f, data = COL.OLD, contrasts = list(fDISCBD =
#> "code_diff", fEW = "contr.treatment"))
#> weights: lw
#> 
#> Moran I statistic standard deviate = 1.5595, p-value = 0.05944
#> alternative hypothesis: greater
#> sample estimates:
#> Observed Moran I      Expectation         Variance 
#>      -0.02000822      -0.23955333       0.01981926 
#> 
if (run) {
# use stats::contr.treatment for both factors
lm.morantest(lm(f, data=COL.OLD, contrasts=list(fDISCBD="contr.treatment", fEW="contr.treatment")),
 lw)
}
#> 
#>  Global Moran I for regression residuals
#> 
#> data:  
#> model: lm(formula = f, data = COL.OLD, contrasts = list(fDISCBD =
#> "contr.treatment", fEW = "contr.treatment"))
#> weights: lw
#> 
#> Moran I statistic standard deviate = 1.3866, p-value = 0.08278
#> alternative hypothesis: greater
#> sample estimates:
#> Observed Moran I      Expectation         Variance 
#>      -0.02000822      -0.55316698       0.14784790 
#> 
if (run) {
# use codingMatrices::code_control for factor
lm.morantest(lm(f, data=COL.OLD, contrasts=list(fDISCBD="code_diff", fEW="code_control")), lw)
}
#> 
#>  Global Moran I for regression residuals
#> 
#> data:  
#> model: lm(formula = f, data = COL.OLD, contrasts = list(fDISCBD =
#> "code_diff", fEW = "code_control"))
#> weights: lw
#> 
#> Moran I statistic standard deviate = 1.568, p-value = 0.05844
#> alternative hypothesis: greater
#> sample estimates:
#> Observed Moran I      Expectation         Variance 
#>      -0.02000822      -0.28007346       0.02750775 
#> 
if (run) {
# use codingMatrices::code_control for both
lm.morantest(lm(f, data=COL.OLD, contrasts=list(fDISCBD="code_control", fEW="code_control")), lw)
}
#> 
#>  Global Moran I for regression residuals
#> 
#> data:  
#> model: lm(formula = f, data = COL.OLD, contrasts = list(fDISCBD =
#> "code_control", fEW = "code_control"))
#> weights: lw
#> 
#> Moran I statistic standard deviate = 1.4316, p-value = 0.07613
#> alternative hypothesis: greater
#> sample estimates:
#> Observed Moran I      Expectation         Variance 
#>      -0.02000822      -0.26418630       0.02909259 
#> 
```
