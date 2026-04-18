# Cluster Classifications for Local Indicators of Spatial Association and Local Indicators for Categorical Data

Used to return a factor showing so-called cluster classification for
local indicators of spatial association for local Moran's I, local
Geary's C (and its multivariate variant) and local Getis-Ord G. This
factor vector can be added to a spatial object for mapping. When `obj`
is of class `licd`, a list of up to six factors for measures of local
composition (analytical and permutation), local configuration
(analytical and permutation), and combined measures, both the
interaction of composition and configuration, and a simplified recoding
of these.

## Usage

``` r
hotspot(obj, ...)

# Default S3 method
hotspot(obj, ...)

# S3 method for class 'localmoran'
hotspot(obj, Prname, cutoff=0.005, quadrant.type="mean",
 p.adjust="fdr", droplevels=TRUE, ...)
# S3 method for class 'summary.localmoransad'
hotspot(obj, Prname, cutoff=0.005,
 quadrant.type="mean", p.adjust="fdr", droplevels=TRUE, ...)
# S3 method for class 'data.frame.localmoranex'
hotspot(obj, Prname, cutoff=0.005,
 quadrant.type="mean", p.adjust="fdr", droplevels=TRUE, ...)

# S3 method for class 'localG'
hotspot(obj, Prname, cutoff=0.005, p.adjust="fdr", droplevels=TRUE, ...)

# S3 method for class 'localC'
hotspot(obj, Prname, cutoff=0.005, p.adjust="fdr", droplevels=TRUE, ...)
# S3 method for class 'licd'
hotspot(obj, type = "both", cutoff = 0.05, p.adjust = "none", 
 droplevels = TRUE, control = list(), ...)
# S3 method for class 'local_jc_uni'
hotspot(obj, cutoff=0.05, p.adjust="none", ...)
```

## Arguments

- obj:

  An object of class `localmoran`, `localC` or `localG`

- Prname:

  A character string, the name of the column containing the probability
  values to be classified by cluster type if found “interesting”

- cutoff:

  Default 0.005, the probability value cutoff larger than which the
  observation is not found “interesting”

- p.adjust:

  Default `"fdr"`, the
  [`p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) method used, one
  of
  `c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")`

- droplevels:

  Default `TRUE`, should empty levels of the input cluster factor be
  dropped

- quadrant.type:

  Default `"mean"`, for `"localmoran"` objects only, can be
  `c("mean", "median", "pysal")` to partition the Moran scatterplot;
  `"mean"` partitions on the means of the variable and its spatial lag,
  `"median"` on medians of the variable and its spatial lag, `"pysal"`
  at zero for the centred variable and its spatial lag

- type:

  When `obj` is of class `licd`, default `both`, may also be `comp` for
  local composition or `config` for local configuration; if `uni`
  replicates local univariate joincount

- control:

  When `obj` is of class `licd`, default `binomial_sidak` 2,
  `binomial_overlap` TRUE, `jcm_sidak` 3. `binomial_overlap` may be set
  FALSE to avoid the Binomial probability values summing to more than
  unity - the tests in Boots (2003, p. 141) do overlap (`>=` and `<=`),
  and the Šidák exponents may be set to 1 to prevent by-observation
  correction for 2 Binomial and 3 Normal probability values per
  observation

- ...:

  other arguments passed to methods.

## Value

A factor showing so-called cluster classification for local indicators
of spatial association. When `obj` is of class `licd`, a list of up to
six factors for measures of local composition (analytical and
permutation), local configuration (analytical and permutation), and
combined measures, both the interaction of composition and
configuration, and a simplified recoding of these.

## See also

[`licd_multi`](https://r-spatial.github.io/spdep/reference/licd_multi.md)

## Author

Roger Bivand

## Examples

``` r
orig <- spData::africa.rook.nb
listw <- nb2listw(orig)
x <- spData::afcon$totcon

set.seed(1)
C <- localC_perm(x, listw)
Ch <- hotspot(C, Prname="Pr(z != E(Ci)) Sim", cutoff=0.05, p.adjust="none")
table(addNA(Ch))
#> 
#> High-High   Low-Low      <NA> 
#>         4         1        37 
set.seed(1)
I <- localmoran_perm(x, listw)
Ih <- hotspot(I, Prname="Pr(z != E(Ii)) Sim", cutoff=0.05, p.adjust="none")
table(addNA(Ih))
#> 
#> High-High      <NA> 
#>         6        36 
Is <- summary(localmoran.sad(lm(x ~ 1), nb=orig))
Ish <- hotspot(Is, Prname="Pr. (Sad)", cutoff=0.05, p.adjust="none")
table(addNA(Ish))
#> 
#> High-High      <NA> 
#>         5        37 
Ie <- as.data.frame(localmoran.exact(lm(x ~ 1), nb=orig))
Ieh <- hotspot(Ie, Prname="Pr. (exact)", cutoff=0.05, p.adjust="none")
table(addNA(Ieh))
#> 
#> High-High      <NA> 
#>         5        37 
set.seed(1)
G <- localG_perm(x, listw)
Gh <- hotspot(G, Prname="Pr(z != E(Gi)) Sim", cutoff=0.05, p.adjust="none")
table(addNA(Gh))
#> 
#> High <NA> 
#>    6   36 
```
