# Rao's score (a.k.a Lagrange Multiplier) diagnostics for spatial dependence in linear models

The function reports the estimates of tests chosen among five statistics
for testing for spatial dependence in linear models. The statistics are
the simple RS test for error dependence (“RSerr”), the simple RS test
for a missing spatially lagged dependent variable (“RSlag”), variants of
these adjusted for the presence of the other (“adjRSerr” tests for error
dependence in the possible presence of a missing lagged dependent
variable, “adjRSlag” the other way round), and a portmanteau test
(“SARMA”, in fact “RSerr” + “adjRSlag”). Note: from spdep 1.3-2, the
tests are re-named “RS” - Rao's score tests, rather than “LM” - Lagrange
multiplier tests to match the naming of tests from the same family in
`SDM.RStests`.

## Usage

    lm.RStests(model, listw, zero.policy=attr(listw, "zero.policy"), test="RSerr",
     spChk=NULL, naSubset=TRUE)
    lm.LMtests(model, listw, zero.policy=attr(listw, "zero.policy"), test="LMerr",
     spChk=NULL, naSubset=TRUE)
    # S3 method for class 'RStestlist'
    print(x, ...)
    # S3 method for class 'RStestlist'
    summary(object, p.adjust.method="none", ...)
    # S3 method for class 'RStestlist.summary'
    print(x, digits=max(3, getOption("digits") - 2), ...)
    <!-- %tracew(listw) -->

## Arguments

- model:

  an object of class `lm` returned by `lm`, or optionally a vector of
  externally calculated residuals (run though `na.omit` if any NAs
  present) for use when only "RSerr" is chosen; weights and offsets
  should not be used in the `lm` object

- listw:

  a `listw` object created for example by `nb2listw`, expected to be
  row-standardised (W-style)

- zero.policy:

  default `attr(listw, "zero.policy")` as set when `listw` was created,
  if attribute not set, use global option value; if TRUE assign zero to
  the lagged value of zones without neighbours, if FALSE assign NA

- test:

  a character vector of tests requested chosen from RSerr, RSlag,
  adjRSerr, adjRSlag, SARMA; test="all" computes all the tests.

- spChk:

  should the data vector names be checked against the spatial objects
  for identity integrity, TRUE, or FALSE, default NULL to use
  [`get.spChkOption()`](https://r-spatial.github.io/spdep/reference/set.spChkOption.md)

- naSubset:

  default TRUE to subset listw object for omitted observations in model
  object (this is a change from earlier behaviour, when the
  `model$na.action` component was ignored, and the listw object had to
  be subsetted by hand)

- x, object:

  object to be printed

- p.adjust.method:

  a character string specifying the probability value adjustment (see
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html)) for multiple
  tests, default "none"

- digits:

  minimum number of significant digits to be used for most numbers

- ...:

  printing arguments to be passed through

## Details

The two types of dependence are for spatial lag \\\rho\\ and spatial
error \\\lambda\\:

\$\$ \mathbf{y} = \mathbf{X \beta} + \rho \mathbf{W\_{(1)} y} +
\mathbf{u}, \$\$ \$\$ \mathbf{u} = \lambda \mathbf{W\_{(2)} u} +
\mathbf{e} \$\$

where \\\mathbf{e}\\ is a well-behaved, uncorrelated error term. Tests
for a missing spatially lagged dependent variable test that \\\rho =
0\\, tests for spatial autocorrelation of the error \\\mathbf{u}\\ test
whether \\\lambda = 0\\. \\\mathbf{W}\\ is a spatial weights matrix; for
the tests used here they are identical.

If factors or ordered factors are among explanatory variables in the
regression, the coding of the contrasts affects the outcomes of
significance tests.

## Value

A list of class `RStestlist` of `htest` objects, each with:

- statistic:

  the value of the Rao's score (a.k.a Lagrange multiplier) test.

- parameter:

  number of degrees of freedom

- p.value:

  the p-value of the test.

- method:

  a character string giving the method used.

- data.name:

  a character string giving the name(s) of the data.

## References

Anselin, L. 1988 Spatial econometrics: methods and models. (Dordrecht:
Kluwer); Anselin, L., Bera, A. K., Florax, R. and Yoon, M. J. 1996
Simple diagnostic tests for spatial dependence. Regional Science and
Urban Economics, 26, 77–104
[doi:10.1016/0166-0462(95)02111-6](https://doi.org/10.1016/0166-0462%2895%2902111-6)
; Malabika Koley (2024) Specification Testing under General Nesting
Spatial Model, <https://sites.google.com/view/malabikakoley/research>.

## Author

Roger Bivand <Roger.Bivand@nhh.no> and Andrew Bernat

## See also

[`lm`](https://rdrr.io/r/stats/lm.html),
[`SD.RStests`](https://r-spatial.github.io/spdep/reference/SD.RStests.md)

## Examples

``` r
data(oldcol)
oldcrime.lm <- lm(CRIME ~ HOVAL + INC, data = COL.OLD)
summary(oldcrime.lm)
#> 
#> Call:
#> lm(formula = CRIME ~ HOVAL + INC, data = COL.OLD)
#> 
#> Residuals:
#>     Min      1Q  Median      3Q     Max 
#> -34.418  -6.388  -1.580   9.052  28.649 
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)  68.6190     4.7355  14.490  < 2e-16 ***
#> HOVAL        -0.2739     0.1032  -2.654   0.0109 *  
#> INC          -1.5973     0.3341  -4.780 1.83e-05 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Residual standard error: 11.43 on 46 degrees of freedom
#> Multiple R-squared:  0.5524, Adjusted R-squared:  0.5329 
#> F-statistic: 28.39 on 2 and 46 DF,  p-value: 9.341e-09
#> 
lw <- nb2listw(COL.nb)
res <- lm.RStests(oldcrime.lm, listw=lw, test="all")
summary(res)
#>  Rao's score (a.k.a Lagrange multiplier) diagnostics for spatial
#>  dependence
#> data:  
#> model: lm(formula = CRIME ~ HOVAL + INC, data = COL.OLD)
#> test weights: lw
#>  
#>          statistic parameter  p.value   
#> RSerr     5.723131         1 0.016743 * 
#> RSlag     9.363684         1 0.002213 **
#> adjRSerr  0.079495         1 0.777983   
#> adjRSlag  3.720048         1 0.053763 . 
#> SARMA     9.443178         2 0.008901 **
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
if (require("spatialreg", quietly=TRUE)) {
  oldcrime.slx <- lmSLX(CRIME ~ HOVAL + INC, data = COL.OLD, listw=lw)
  summary(lm.RStests(oldcrime.slx, listw=lw, test=c("adjRSerr", "adjRSlag")))
}
#>  Rao's score (a.k.a Lagrange multiplier) diagnostics for spatial
#>  dependence
#> data:  
#> model: lm(CRIME ~ HOVAL + INC + lag.HOVAL + lag.INC, data = COL.OLD,
#> listw = lw)
#> test weights: lw
#>  
#>              statistic parameter p.value
#> GNM_adjRSerr   0.10933         1  0.7409
#> GNM_adjRSlag   0.61580         1  0.4326
run <- require("codingMatrices", quietly=TRUE)
if (run) {
COL.OLD$fEW <- factor(COL.OLD$EW)
COL.OLD$fDISCBD <- ordered(cut(COL.OLD$DISCBD, c(0, 1.5, 3, 4.5, 6)))
f <- formula(CRIME ~ INC + HOVAL + fDISCBD*fEW)
lw <- nb2listw(COL.nb, style="W")
# default codings
summary(lm.RStests(lm(f, data=COL.OLD, contrasts=list(fDISCBD="contr.poly", fEW="contr.treatment")),
 lw, test=c("adjRSerr", "adjRSlag")))
}
#>  Rao's score (a.k.a Lagrange multiplier) diagnostics for spatial
#>  dependence
#> data:  
#> model: lm(formula = f, data = COL.OLD, contrasts = list(fDISCBD =
#> "contr.poly", fEW = "contr.treatment"))
#> test weights: lw
#>  
#>          statistic parameter p.value
#> adjRSerr   0.98727         1  0.3204
#> adjRSlag   1.14732         1  0.2841
if (run) {
# use codingMatrices::code_diff for ordered factor
summary(lm.RStests(lm(f, data=COL.OLD, contrasts=list(fDISCBD="code_diff", fEW="contr.treatment")),
 lw, test=c("adjRSerr", "adjRSlag")))
}
#>  Rao's score (a.k.a Lagrange multiplier) diagnostics for spatial
#>  dependence
#> data:  
#> model: lm(formula = f, data = COL.OLD, contrasts = list(fDISCBD =
#> "code_diff", fEW = "contr.treatment"))
#> test weights: lw
#>  
#>          statistic parameter p.value
#> adjRSerr   0.03415         1  0.8534
#> adjRSlag  -0.01572         1  1.0000
if (run) {
# use codingMatrices::code_control for both
summary(lm.RStests(lm(f, data=COL.OLD, contrasts=list(fDISCBD="code_control", fEW="code_control")),
 lw, test=c("adjRSerr", "adjRSlag")))
}
#>  Rao's score (a.k.a Lagrange multiplier) diagnostics for spatial
#>  dependence
#> data:  
#> model: lm(formula = f, data = COL.OLD, contrasts = list(fDISCBD =
#> "code_control", fEW = "code_control"))
#> test weights: lw
#>  
#>           statistic parameter p.value
#> adjRSerr  0.0383759         1  0.8447
#> adjRSlag -0.0062517         1  1.0000
```
