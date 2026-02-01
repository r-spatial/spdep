# Local Indicators for Categorical Data

Local indicators for categorical data combine a measure of local
composition in a window given by the per-observation set of neighbouring
observations, with a local multi-category joincount test simplified to
neighbours with the same or different categories compared to the focal
observation

## Usage

``` r
licd_multi(fx, listw, zero.policy = attr(listw, "zero.policy"), adjust.n = TRUE,
 nsim = 0L, iseed = NULL, no_repeat_in_row = FALSE, control = list())
```

## Arguments

- fx:

  a factor with two or more categories, of the same length as the
  neighbours and weights objects in listw; use of an ordered factor is
  not well understood

- listw:

  a `listw` object created for example by `nb2listw`

- zero.policy:

  default `attr(listw, "zero.policy")` as set when `listw` was created,
  if attribute not set, use global option value; if TRUE assign zero to
  the lagged value of zones without neighbours, if FALSE assign NA

- adjust.n:

  default TRUE, if FALSE the number of observations is not adjusted for
  no-neighbour observations, if TRUE, the number of observations is
  adjusted

- nsim:

  default 0, number of conditonal permutation simulations

- iseed:

  default NULL, used to set the seed; the output will only be
  reproducible if the count of CPU cores across which computation is
  distributed is the same

- no_repeat_in_row:

  default `FALSE`, if `TRUE`, sample conditionally in each row without
  replacements to avoid duplicate values,
  <https://github.com/r-spatial/spdep/issues/124>

- control:

  comp_binary=`TRUE` default TRUE, ignoring other weights styles than
  binary for composition measures,
  binomial_punif_alternative=`"greater"`,
  jcm_same_punif_alternative=`"less"`,
  jcm_diff_punif_alternative=`"greater"`,
  uni_jc_same_punif_alternative=`"two.sided"`, rank_ties.method=`"min"`
  default "min", na.last=`"keep"` leading to rank NA being returned if
  the observed joincount variance is non-positive; if TRUE joincount NAs
  are ranked highest when using rank, for others see ?rank,
  check_reps=`FALSE`, unique_ceiling=1/3 used if check_reps TRUE,
  check_reps=FALSE if TRUE, check and report how many unique draws are
  made among the nsim draws, and if the number of unique draws is less
  than unique_ceiling, compute measures only for unique draws and copy
  out to replicated draws, pysal_rank=`FALSE` use rank with
  rank_ties.method and na.last; if TRUE, use pysal-style ranking to find
  the rank of sum(sims \<= obs, na.rm=TRUE)+1 for pysal_sim_obs "GE",
  sum(sims \< obs, na.rm=TRUE)+1 for the count of observed greater than
  "GT" the simulated values, pysal_sim_obs="GT" may also be "GE",
  xtras=`FALSE` if TRUE return calculated compostion values of BW
  chi-squared, k-colour chi-squared, BW Anscombe, and emulates local
  univariate joincount (requires conditional permutation)

## Details

The original code may be found at
[doi:10.5281/zenodo.4283766](https://doi.org/10.5281/zenodo.4283766)

## Value

- local_comp:

  data.frame object with LICD local composition columns: ID, category_i,
  count_like_i, prop_i, count_nbs_i, pbinom_like_BW, pbinom_unlike_BW,
  pbinom_unlike_BW_alt, chi_BW_i, chi_K_i, anscombe_BW

- local_config:

  data.frame object with LICD local configuration columns: ID,
  jcm_chi_obs, jcm_count_BB_obs, jcm_count_BW_obs, jcm_count_WW_obs,
  pval_jcm_obs_BB, pval_jcm_obs_WW, pval_jcm_obs_BW, only_i_jc

- local_comp_sim:

  data.frame object with permutation-based LICD local composition
  columns: ID, pbinom_like_BW, pbinom_unlike_BW, pbinom_unlike_BW_alt,
  rank_sim_chi_BW, rank_sim_chi_K, rank_sim_anscombe

- local_config_sim:

  data.frame object with permutation-based LICD local configuration
  columns: ID, jcm_chi_sim_rank, pval_jcm_obs_BB, pval_jcm_obs_BW,
  pval_jcm_obs_WW

## References

Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 20;

Upton, G., Fingleton, B. 1985 Spatial data analysis by example: point
pattern and qualitative data, Wiley, pp. 158–170;

Boots, B., 2003. Developing local measures of spatial association for
categorical data. Journal of Geographical Systems 5, 139–160;

Boots, B., 2006. Local configuration measures for categorical spatial
data: binary regular lattices. Journal of Geographical Systems 8 (1),
1–24;

Pietrzak, M.B., Wilk, J., Kossowski, T., Bivand, R.S., 2014. The
application of local indicators for categorical data (LICD) in the
spatial analysis of economic development. Comparative Economic Research
17 (4), 203–220
[doi:10.2478/cer-2014-0041](https://doi.org/10.2478/cer-2014-0041) ;

Bivand, R.S., Wilk, J., Kossowski, T., 2017. Spatial association of
population pyramids across Europe: The application of symbolic data,
cluster analysis and join-count tests. Spatial Statistics 21 (B),
339–361
[doi:10.1016/j.spasta.2017.03.003](https://doi.org/10.1016/j.spasta.2017.03.003)
;

Francesco Carrer, Tomasz M. Kossowski, Justyna Wilk, Michał B. Pietrzak,
Roger S. Bivand, The application of Local Indicators for Categorical
Data (LICD) to explore spatial dependence in archaeological spaces,
Journal of Archaeological Science, 126, 2021,
[doi:10.1016/j.jas.2020.105306](https://doi.org/10.1016/j.jas.2020.105306)

## Author

Roger Bivand <Roger.Bivand@nhh.no> based on earlier code by Tomasz M.
Kossowski, Justyna Wilk and Michał B. Pietrzak

## Note

In order to increase the numbers of neighbours using
[`nblag`](https://r-spatial.github.io/spdep/reference/nblag.md) and
[`nblag_cumul`](https://r-spatial.github.io/spdep/reference/nblag.md) is
advisable; use of binary weights is advised and are in any case used for
the composition measure

## See also

[`joincount.multi`](https://r-spatial.github.io/spdep/reference/joincount.multi.md)

## Examples

``` r
columbus <- st_read(system.file("shapes/columbus.gpkg", package="spData")[1], quiet=TRUE)
HICRIME <- cut(columbus$CRIME, breaks=c(0,35,80), labels=c("low","high"))
(nb <- poly2nb(columbus))
#> Neighbour list object:
#> Number of regions: 49 
#> Number of nonzero links: 236 
#> Percentage nonzero weights: 9.829238 
#> Average number of links: 4.816327 
lw <- nb2listw(nblag_cumul(nblag(nb, 2)), style="B")
obj <- licd_multi(HICRIME, lw)
str(obj)
#> List of 5
#>  $ local_comp      :'data.frame':    49 obs. of  11 variables:
#>   ..$ ID                  : int [1:49] 1 2 3 4 5 6 7 8 9 10 ...
#>   ..$ category_i          : num [1:49] 1 1 1 1 2 1 1 2 1 1 ...
#>   ..$ count_like_i        : num [1:49] 4 4 6 7 12 7 2 8 10 8 ...
#>   ..$ prop_i              : num [1:49] 0.51 0.51 0.51 0.51 0.49 ...
#>   ..$ count_nbs_i         : num [1:49] 5 6 11 14 22 14 11 14 25 17 ...
#>   ..$ pbinom_like_BW      : num [1:49] 0.965 0.881 0.702 0.575 0.769 ...
#>   ..$ pbinom_unlike_BW    : num [1:49] 0.0346 0.1192 0.2979 0.4255 0.2313 ...
#>   ..$ pbinom_unlike_BW_alt: num [1:49] 0.201 0.363 0.528 0.634 0.379 ...
#>   ..$ chi_BW_i            : num [1:49] NA NA NA NA NA NA NA NA NA NA ...
#>   ..$ chi_K_i             : num [1:49] NA NA NA NA NA NA NA NA NA NA ...
#>   ..$ anscombe_BW         : num [1:49] NA NA NA NA NA NA NA NA NA NA ...
#>  $ local_config    :'data.frame':    49 obs. of  9 variables:
#>   ..$ ID              : int [1:49] 1 2 3 4 5 6 7 8 9 10 ...
#>   ..$ jcm_chi_obs     : num [1:49] 0 0.0625 0.625 2.7468 14.1889 ...
#>   ..$ jcm_count_BB_obs: num [1:49] 6 6 11 12 55 15 1 27 30 25 ...
#>   ..$ jcm_count_BW_obs: num [1:49] 4 7 23 31 50 36 15 31 73 53 ...
#>   ..$ jcm_count_WW_obs: num [1:49] 0 1 10 20 21 19 32 9 61 24 ...
#>   ..$ pval_jcm_obs_BB : num [1:49] 1 0.207108 0.70232 0.820405 0.000214 ...
#>   ..$ pval_jcm_obs_WW : num [1:49] 1 0.3946 0.0981 0.0243 0.7839 ...
#>   ..$ pval_jcm_obs_BW : num [1:49] 1.00 1.75e-01 2.08e-01 5.10e-02 3.52e-06 ...
#>   ..$ only_i_jc       : num [1:49] 3 3 5 6 11 6 1 7 9 7 ...
#>  $ local_comp_sim  : NULL
#>  $ local_config_sim: NULL
#>  $ local_uni_sim   : NULL
#>  - attr(*, "timings")=List of 3
#>   ..$ set_up        : 'proc_time' Named num [1:5] 0.001 0 0.001 0 0
#>   .. ..- attr(*, "names")= chr [1:5] "user.self" "sys.self" "elapsed" "user.child" ...
#>   ..$ processing    : 'proc_time' Named num [1:5] 0.092 0 0.092 0 0
#>   .. ..- attr(*, "names")= chr [1:5] "user.self" "sys.self" "elapsed" "user.child" ...
#>   ..$ postprocessing: 'proc_time' Named num [1:5] 0.001 0 0.001 0 0
#>   .. ..- attr(*, "names")= chr [1:5] "user.self" "sys.self" "elapsed" "user.child" ...
#>  - attr(*, "out")= num [1:49, 1:37] 1 1 1 1 2 1 1 2 1 1 ...
#>   ..- attr(*, "ncpus")= int 1
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:37] "category_i" "count_like_i" "prop_i" "count_nbs_i" ...
#>  - attr(*, "ncpus")= int 1
#>  - attr(*, "nsim")= int 0
#>  - attr(*, "con")=List of 12
#>   ..$ comp_binary                  : logi TRUE
#>   ..$ binomial_punif_alternative   : chr "greater"
#>   ..$ jcm_same_punif_alternative   : chr "less"
#>   ..$ jcm_diff_punif_alternative   : chr "greater"
#>   ..$ uni_jc_same_punif_alternative: chr "two.sided"
#>   ..$ rank_ties.method             : chr "min"
#>   ..$ unique_ceiling               : num 0.333
#>   ..$ check_reps                   : logi FALSE
#>   ..$ pysal_rank                   : logi FALSE
#>   ..$ pysal_sim_obs                : chr "GT"
#>   ..$ na.last                      : chr "keep"
#>   ..$ xtras                        : logi FALSE
#>  - attr(*, "uni_jcs")= num [1:49, 1:2] 3 3 5 6 0 6 1 0 9 7 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:49] "1" "2" "3" "4" ...
#>   .. ..$ : chr [1:2] "fxlow" "fxhigh"
#>  - attr(*, "class")= chr [1:2] "licd" "list"
h_obj <- hotspot(obj)
str(h_obj)
#> List of 10
#>  $ ID              : int [1:49] 1 2 3 4 5 6 7 8 9 10 ...
#>  $ local_comp      : Factor w/ 2 levels "Cluster","Dispersed": 2 2 2 2 2 2 2 2 2 2 ...
#>  $ local_comp_sim  : NULL
#>  $ local_config    : Factor w/ 3 levels "Cluster","Dispersed",..: 3 3 3 3 2 3 3 2 3 3 ...
#>  $ local_config_sim: NULL
#>  $ local_uni_sim   : NULL
#>  $ both            : Factor w/ 6 levels "Cluster.Cluster",..: 6 6 6 6 4 6 6 4 6 6 ...
#>  $ both_sim        : NULL
#>  $ both_recode     : Factor w/ 4 levels "Clump","Cluster",..: 4 4 4 4 3 4 4 3 4 4 ...
#>  $ both_recode_sim : NULL
table(h_obj$both_recode)
#> 
#>      Clump    Cluster  Dispersed No cluster 
#>          8          3         11         27 
columbus$both <- h_obj$both_recode
plot(columbus[, "both"])

GDAL37 <- numeric_version(unname(sf::sf_extSoftVersion()["GDAL"]), strict=FALSE)
(GDAL37 <- ifelse(is.na(GDAL37), FALSE, GDAL37 >= "3.7.0"))
#> [1] TRUE
file <- "etc/shapes/GB_2024_southcoast_50m.gpkg.zip"
zipfile <- system.file(file, package="spdep")
if (GDAL37) {
    sc50m <- st_read(zipfile)
} else {
    td <- tempdir()
    bn <- sub(".zip", "", basename(file), fixed=TRUE)
    target <- unzip(zipfile, files=bn, exdir=td)
    sc50m <- st_read(target)
}
#> Reading layer `GB_2024_southcoast_50m' from data source 
#>   `/tmp/RtmpEjBomy/temp_libpath486f623f1b6ee/spdep/etc/shapes/GB_2024_southcoast_50m.gpkg.zip' 
#>   using driver `GPKG'
#> Simple feature collection with 119 features and 19 fields
#> Geometry type: MULTIPOLYGON
#> Dimension:     XY
#> Bounding box:  xmin: 82643.12 ymin: 5342.9 xmax: 640301.6 ymax: 187226.2
#> Projected CRS: OSGB36 / British National Grid
sc50m$Winner <- factor(sc50m$Winner, levels=c("Con", "Green", "Lab", "LD"))
plot(sc50m[,"Winner"], pal=c("#2297E6", "#61D04F", "#DF536B", "#F5C710"))

nb_sc_50m <- poly2nb(sc50m, row.names=as.character(sc50m$Constituency))
#> Warning: neighbour object has 2 sub-graphs;
#> if this sub-graph count seems unexpected, try increasing the snap argument.
sub2 <- attr(nb_sc_50m, "region.id")[attr(nb_sc_50m, "ncomp")$comp.id == 2L]
iowe <- match(sub2[1], attr(nb_sc_50m, "region.id"))
diowe <- c(st_distance(sc50m[iowe,], sc50m))
meet_criterion <- sum(diowe <= units::set_units(5000, "m"))
cands <- attr(nb_sc_50m, "region.id")[order(diowe)[1:meet_criterion]]
nb_sc_50m_iowe <- addlinks1(nb_sc_50m, from = cands[1],
 to = cands[3:meet_criterion])
ioww <- match(sub2[2], attr(nb_sc_50m, "region.id"))
dioww <- c(st_distance(sc50m[ioww,], sc50m))
meet_criterion <- sum(dioww <= units::set_units(5000, "m"))
cands <- attr(nb_sc_50m, "region.id")[order(dioww)[1:meet_criterion]]
nb_sc_50m_iow <- addlinks1(nb_sc_50m_iowe, from = cands[2], to = cands[3:meet_criterion])
nb_sc_1_2 <- nblag_cumul(nblag(nb_sc_50m_iow, 2))
lw <- nb2listw(nb_sc_1_2, style="B")
licd_obj <- licd_multi(sc50m$Winner, lw)
h_obj <- hotspot(licd_obj)
sc50m$both <- h_obj$both_recode
plot(sc50m[, "both"])

ljc <- local_joincount_uni(factor(sc50m$Winner == "LD"), chosen="TRUE", lw)
sc50m$LD_pv <- ljc[, 2]
plot(sc50m[, "LD_pv"], breaks=c(0, 0.025, 0.05, 0.1, 0.5, 1))
```
