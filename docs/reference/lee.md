# Compute Lee's statistic

A simple function to compute Lee's L statistic for bivariate spatial
data; \$\$L(x,y) = \frac{n}{\sum\_{i=1}^{n}(\sum\_{j=1}^{n}w\_{ij})^2}
\frac{\sum\_{i=1}^{n}(\sum\_{j=1}^{n}w\_{ij}(x_i-\bar{x}))
((\sum\_{j=1}^{n}w\_{ij}(y_j-\bar{y}))}{\sqrt{\sum\_{i=1}^{n}(x_i -
\bar{x})^2} \sqrt{\sum\_{i=1}^{n}(y_i - \bar{y})^2}} \$\$

## Usage

``` r
lee(x, y, listw, n, S2, zero.policy=attr(listw, "zero.policy"), NAOK=FALSE)
```

## Arguments

- x:

  a numeric vector the same length as the neighbours list in listw

- y:

  a numeric vector the same length as the neighbours list in listw

- listw:

  a `listw` object created for example by `nb2listw`

- n:

  number of zones

- S2:

  Sum of squared sum of weights by rows.

- zero.policy:

  default `attr(listw, "zero.policy")` as set when `listw` was created,
  if attribute not set, use global option value; if TRUE assign zero to
  the lagged value of zones without neighbours, if FALSE assign NA

- NAOK:

  if 'TRUE' then any 'NA' or 'NaN' or 'Inf' values in x are passed on to
  the foreign function. If 'FALSE', the presence of 'NA' or 'NaN' or
  'Inf' values is regarded as an error.

## Value

a list of

- L:

  Lee's L statistic

- local L:

  Lee's local L statistic

## References

Lee (2001). Developing a bivariate spatial association measure: An
integration of Pearson's r and Moran's I. J Geograph Syst 3: 369-385

## Author

Roger Bivand and Virgiio GÃ³mez-Rubio <Virgilio.Gomez@uclm.es>

## See also

[`lee.mc`](https://r-spatial.github.io/spdep/reference/lee.mc.md)

## Examples

``` r
data(boston, package="spData")
lw<-nb2listw(boston.soi)

x<-boston.c$CMEDV
y<-boston.c$CRIM
z<-boston.c$RAD

Lxy<-lee(x, y, lw, length(x), zero.policy=TRUE)
Lxz<-lee(x, z, lw, length(x), zero.policy=TRUE)
```
