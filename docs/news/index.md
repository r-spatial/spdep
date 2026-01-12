# Changelog

## Version 1.4-2 (development)

## Version 1.4-1 (2025-08-31)

CRAN release: 2025-08-31

- add data argument to `SD.RStests` for formula Durbin case

## Version 1-3-13 (2025-06-10)

- Adds Bavaud’s multivariate `spatialdelta` with support functions and
  methods

- Add note on changes to output from tests for error autocorrelation if
  contrast codings are set to non-default values

## Version 1.3-11 (2025-04-24)

CRAN release: 2025-04-24

- introduce warnings for factors (categorical variables) in Durbin
  models (`SD.RStests`)

- remove \|\> in vignette to avoid R \>= 4.1 dependency
  <https://stat.ethz.ch/pipermail/r-devel/2025-January/083768.html>

- update reference to Koley (2024) in man/SD.RStests.Rd

- add warnings to Durbin terms including categorical variables in
  `SD.RStests`

## Version 1.3-10 (2025-01-20)

CRAN release: 2025-01-20

- `R_NO_REMAP` problem with R \< 4.4.1 fixed by replacing
  `COPY_TO_USER_STRING` by `Rf_mkChar`, because include/Rdefines.h
  before revision 86416 set `COPY_TO_USER_STRING` as `mkChar` which is
  not defined as `Rf_mkChar` when `R_NO_REMAP` is defined,
  [\#176](https://github.com/r-spatial/spdep/issues/176), thanks to
  Edzer Pebesma

- disambiguate which `skater` in examples to satisfy `pkgdown`, which
  used
  [`rgeoda::skater`](https://geodacenter.github.io/rgeoda/reference/skater.html);
  other patches to examples to satisfy `pkgdown`

- use `inherits` in `skater`

## Version 1.3-9 (2025-01-16)

CRAN release: 2025-01-16

- revisit `diffnb` and set operations like `union.nb` and `setdiff.nb`
  following up [\#175](https://github.com/r-spatial/spdep/issues/175);
  `diffnb` largely rewritten and should no longer generate deformed
  output; set operations modified to match base functions actions

- convert `error` to `Rf_error`, `length` to `Rf_length` etc. to
  accommodate `R_NO_REMAP`, see
  <https://github.com/r-spatial/spdep/commit/b61f6b17be09383c94b45e912f9213735aa62212>
  for R 4.5

- re-instate **rgeoda** references

## Version 1.3-8 (2024-12-02)

CRAN release: 2024-12-02

- Remove remaining `spData` ESRI shapefile use

## Version 1.3-7 (2024-11-25)

CRAN release: 2024-11-25

- (temporarily) remove **rgeoda** references until it is successfully
  re-submitted to CRAN

- add `write.swmdbf`
  [\#171](https://github.com/r-spatial/spdep/issues/171) to complement
  [\#163](https://github.com/r-spatial/spdep/issues/163)

- modify defaults for `licd_multi`

## Version 1.3-6 (2024-09-13)

CRAN release: 2024-09-13

- adding vignette describing recent changes in `poly2nb` from
  [\#162](https://github.com/r-spatial/spdep/issues/162), subgraph and
  no-neighbour (island) handling

- adding prototype of LICD ESDA function `licd_multi` and `hotspot`
  method

- add `read.swmdbf2listw`
  [\#163](https://github.com/r-spatial/spdep/issues/163) for reading DBF
  files exported from ArcGIS representing SWM objects; note that there
  will be problems when the observation IDs are not known, see help file

- [\#162](https://github.com/r-spatial/spdep/issues/162) add option for
  no-neighbour checking for `poly2nb` - default report whether
  no-neighbour observations are present

- [\#162](https://github.com/r-spatial/spdep/issues/162) change the
  default `snap=` argument to `poly2nb` to 10mm

- Condition on forthcoming `tmap` 4

- [\#160](https://github.com/r-spatial/spdep/issues/160) handle
  `n.comp.nb` delay in `print.nb` and elsewhere when the total number of
  neighbours is large

## Version 1.3-5 (2024-06-10)

CRAN release: 2024-06-10

- [\#157](https://github.com/r-spatial/spdep/issues/157) migrate ESRI
  Shapefile to GPKG files; convert bhicv.shp to GPKG

- [\#155](https://github.com/r-spatial/spdep/issues/155) Throw error if
  `hotspot` despatched on object without a `"quadr"` attribute

- [\#154](https://github.com/r-spatial/spdep/issues/154) turn on `s2` in
  vignette

## Version 1.3-4 (2024-05-31)

CRAN release: 2024-05-31

- add `scale` argument to `geary.test`, `geary.mc` and `geary`
  [\#151](https://github.com/r-spatial/spdep/issues/151), and
  appropriate tests

- Introduce error in `knearneigh` for `k` less than the largest count of
  identical points; if encountered, increase `k`

- remove spurious warning in `knearneigh` for longlat geometries

- fix <https://github.com/edzer/sdsr/issues/121>, wrong assignment of
  old test names in `lmRStests`

- fix [\#144](https://github.com/r-spatial/spdep/issues/144) in
  `plot.nb` and `nb2lines`

## Version 1.3-3 (2024-02-07)

CRAN release: 2024-02-07

- change `lm.LMtests` to `lm.RStests` and re-name Lagrange multiplier to
  Rao’s score; add `GNM_` prefix to test names if the input object
  inherits from `SlX` created by
  [`spatialreg::lmSLX`](https://r-spatial.github.io/spatialreg/reference/SLX.html)
  (Koley, forthcoming)

- add `SD.RStests` implementation of Rao’s score tests for spatial
  Durbin models (Koley and Bera, 2024) and for SDEM models (Koley,
  forthcoming)

- [\#143](https://github.com/r-spatial/spdep/issues/143) `row.names`
  pass-through in `poly2nb` corrected, harmonised `row.names`
  pass-through also in `nbdists` and `dnearneigh`

- [\#139](https://github.com/r-spatial/spdep/issues/139) add `na.action`
  argument to `geary.test`, `geary.mc` and `globalG.test`

- add `style` to `sn2listw` use in `tri2nb`

## Version 1.3-1 (2023-11-23)

CRAN release: 2023-11-23

- functions creating `nb` objects now warn if the object has a sub-graph
  count of \> 1 and `get.SubgraphOption` is `TRUE` (default `FALSE`):
  `complement.nb`, `diffnb`, `dnearneigh`, `droplinks`, `edit.nb`,
  `graph2nb`, `knn2nb`, `nb2blocknb`, `nblag`, `nblag_cumul`, `poly2nb`,
  `read.gal`, `read.gwt2nb`, `setdiff.nb`, `tolerance.nb`, `tri2nb`,
  `union.nb`

- `summary.nb`, `print.nb`, `summary.listw` and `print.listw` now report
  the subgraph count from `n.comp.nb` if it is more than one

- `subset.nb` now reports if the subgraph count of the neighbour object
  increases on subsetting

- adding a `zero.policy` attribute to functions creating `listw`
  objects: `nb2listw`, `sn2listw`, `mat2listw`, `nb2listwdist`. Default
  `zero.policy=` argument updated to use `attr(., "zero.policy")` in
  `summary.listw`, `print.listw`, `moran`, `moran.test`, `moran.mc`,
  `moran.plot`, `geary.mc`, `geary`, `geary.test`, `globalG.test`,
  `joincount.test`, `joincount.mc`, `joincount.multi`, `localC`,
  `localC_perm`, `localmoran`, `localmoran_perm`, `localG`,
  `localG_perm`, `lee`, `lee.test`, `lee.mc`, `lm.morantest`,
  `lm.LMtests`, `sp.mantel.mc`, `listw2star`, `lag.listw`,
  `lm.morantest`, `lm.LMtests`, `subset.listw`, `EBImoran.mc`, `LOSH`,
  `LOSH.mc`, `LOSH.cs`, `lm.morantest.exact`and `lm.morantest.sad`

- confusing error message in
  [`moran.plot()`](https://r-spatial.github.io/spdep/reference/moran.plot.md)
  if no-neighbour cases, but `zero.policy=FALSE`

- replace `rgrass7` with `rgrass` in vignette

- fix [\#133](https://github.com/r-spatial/spdep/issues/133) (`edit.nb`
  affected by not attaching `sp`)

## Version 1.2-8 (2023-02-28)

CRAN release: 2023-02-28

- [`mat2listw()`](https://r-spatial.github.io/spdep/reference/mat2listw.md)
  warning if no `style=` argument given, or if `M"` is given
  <https://github.com/r-spatial/spatialreg/issues/24>,
  <https://github.com/r-spatial/spatialreg/issues/23>.

- remaining users of `run_perm()` -
  [`localC()`](https://r-spatial.github.io/spdep/reference/localC.md),
  [`localmoran_bv()`](https://r-spatial.github.io/spdep/reference/localmoran_bv.md)
  and
  [`local_joincount_uni()`](https://r-spatial.github.io/spdep/reference/local_joincount_uni.md)
  get `no_repeat_in_row=` arguments.

- Address 2) in [\#124](https://github.com/r-spatial/spdep/issues/124);
  [`localG_perm()`](https://r-spatial.github.io/spdep/reference/localG.md)
  and
  [`localmoran_perm()`](https://r-spatial.github.io/spdep/reference/localmoran.md)
  get `no_repeat_in_row=` arguments to use conditional permutation
  without replacement by sample vectors; the default implementation uses
  sampling with replacement, which is acceptable across simulation
  draws, but arguably less acceptable within draws. Feedback would be
  valued.

- Address 1) in [\#124](https://github.com/r-spatial/spdep/issues/124);
  [`localG_perm()`](https://r-spatial.github.io/spdep/reference/localG.md)
  and
  [`localG()`](https://r-spatial.github.io/spdep/reference/localG.md)
  now return the same analytical standard deviates. The standard
  deviates from the simulated distributions are now returned in
  `attr(., "internals")[,"StdDev.Gi"]` from
  [`localG_perm()`](https://r-spatial.github.io/spdep/reference/localG.md),
  as are p-values, etc.

- move **sp** from Depends to Imports, to reduce the visual impression
  that **sp** objects are required for **spdep**; **sf** objects are now
  preferred, but **sp** objects can be used as before, although users
  may need to attach **sp** expliciitly.

- fix [\#121](https://github.com/r-spatial/spdep/issues/121) and
  [\#123](https://github.com/r-spatial/spdep/issues/123); correcting
  returned values for
  [`localG_perm()`](https://r-spatial.github.io/spdep/reference/localG.md)
  when estimating the G-star measure (fix self x values and weights)

- address [\#120](https://github.com/r-spatial/spdep/issues/120), moving
  documentation of
  [`listw2U()`](https://r-spatial.github.io/spdep/reference/nb2listw.md)
  from
  [`?lm.morantest`](https://r-spatial.github.io/spdep/reference/lm.morantest.md)
  to
  [`?nb2listw`](https://r-spatial.github.io/spdep/reference/nb2listw.md)

- addressing [\#119](https://github.com/r-spatial/spdep/issues/119) for
  interpretation of
  [`moran_bv()`](https://r-spatial.github.io/spdep/reference/moran_bv.md)
  results

- adding [\#116](https://github.com/r-spatial/spdep/issues/116), René
  Westerholt

- fix [\#113](https://github.com/r-spatial/spdep/issues/113), too low R
  version for
  [`grDevices::hcl.colors()`](https://rdrr.io/r/grDevices/palettes.html)

- addressing [\#111](https://github.com/r-spatial/spdep/issues/111) by
  Josiah Parry

- PRs from René Westerholt, ending with
  [\#109](https://github.com/r-spatial/spdep/issues/109), for the local
  GS measure

## Version 1.2-7 (2022-10-01)

CRAN release: 2022-10-01

- [\#103](https://github.com/r-spatial/spdep/issues/103) refactoring
  local joincount test by Josiah Parry

- add `hotspot` methods for `localmoran` (analytical, permutation,
  Saddlepoint and exact), `localC` (univariate and multivariate) and
  `localG` (analytical and permutation)

- [\#95](https://github.com/r-spatial/spdep/issues/95) add `"two.sided"`
  to
  [`lee.mc()`](https://r-spatial.github.io/spdep/reference/lee.mc.md),
  same for
  [`sp.mantel.mc()`](https://r-spatial.github.io/spdep/reference/sp.mantel.mc.md),
  [`EBImoran.mc()`](https://r-spatial.github.io/spdep/reference/EBImoran.mc.md),
  [`joincount.mc()`](https://r-spatial.github.io/spdep/reference/joincount.mc.md),
  [`geary.mc()`](https://r-spatial.github.io/spdep/reference/geary.mc.md)

- [\#92](https://github.com/r-spatial/spdep/issues/92),
  [\#93](https://github.com/r-spatial/spdep/issues/93),
  [\#94](https://github.com/r-spatial/spdep/issues/94),
  [\#96](https://github.com/r-spatial/spdep/issues/96),
  [\#97](https://github.com/r-spatial/spdep/issues/97) contributions of
  prototype bivariate Moran, local bivariate Moran and local joincount
  and bivariate joincount tests by Josiah Parry

- [\#91](https://github.com/r-spatial/spdep/issues/91)
  [`tolerance.nb()`](https://r-spatial.github.io/spdep/reference/tolerance.nb.md)
  update by F. Guillaume Blanchet

- updating coercion for **Matrix** 1.4-2

- fix ncpus issue in dontrun examples

- remove suggested packages **rgdal**, **rgeos**, **maptools**

## Version 1.2-5 (2022-08-11)

CRAN release: 2022-08-11

- permit use of data.frame or tibble as matrix for functions creating
  neighbour objects from 2D points (preferred use an object inheriting
  from \``"SpatialPoints"` or `"sfc"`)

- fix [\#87](https://github.com/r-spatial/spdep/issues/87) wrong logic
  in infinite weights in
  [`nb2listwdist()`](https://r-spatial.github.io/spdep/reference/nb2listwdist.md)

- <https://github.com/r-spatial/s2/pull/174> speeds up
  [`dnearneigh()`](https://r-spatial.github.io/spdep/reference/dnearneigh.md)
  for geographical coordinates without using
  [`s2::s2_closest_edges()`](https://r-spatial.github.io/s2/reference/s2_closest_feature.html).

- Adapting vignettes for absence of **rgdal** and **maptools**

## Version 1.2-4 (2022-04-18)

CRAN release: 2022-04-18

- added
  [`remove.self()`](https://r-spatial.github.io/spdep/reference/include.self.md),
  thanks to Josiah Parry
  [\#83](https://github.com/r-spatial/spdep/issues/83).

- unescape underscores in help pages.

## Version 1.2-3 (2022-03-29)

CRAN release: 2022-03-29

- replace deprecated S-compatibility macros `DOUBLE_`

- [\#81](https://github.com/r-spatial/spdep/issues/81) improved
  [`dnearneigh()`](https://r-spatial.github.io/spdep/reference/dnearneigh.md)
  help page.

- [\#79](https://github.com/r-spatial/spdep/issues/79) remove `"htest"`
  class from
  [`LOSH.mc()`](https://r-spatial.github.io/spdep/reference/LOSH.mc.md)
  output object.

- Added GA SI article to citations.

## Version 1.2-2 (2022-01-28)

CRAN release: 2022-01-28

- Replace [`rainbow()`](https://rdrr.io/r/grDevices/palettes.html) by
  `hcl.colors(..., "Set 2")` in
  [`plot.skater()`](https://r-spatial.github.io/spdep/reference/plot.skater.md).

- Add link to R-sig-geo thread on
  [`EBlocal()`](https://r-spatial.github.io/spdep/reference/EBlocal.md)
  NaN estimates when many counts are zero on help page.

- Revise and add documentation for object returned by
  [`localC_perm()`](https://r-spatial.github.io/spdep/reference/localC.md)
  [\#68](https://github.com/r-spatial/spdep/issues/68)
  [\#72](https://github.com/r-spatial/spdep/issues/72)
  [\#73](https://github.com/r-spatial/spdep/issues/73)
  [\#74](https://github.com/r-spatial/spdep/issues/74)
  [\#75](https://github.com/r-spatial/spdep/issues/75)
  [\#76](https://github.com/r-spatial/spdep/issues/76)

- [`localmoran.sad()`](https://r-spatial.github.io/spdep/reference/localmoran.sad.md),
  [`localmoran.exact()`](https://r-spatial.github.io/spdep/reference/localmoran.exact.md)
  and
  [`localmoran.exact.alt()`](https://r-spatial.github.io/spdep/reference/localmoran.exact.md)
  will now use multiple compute nodes if needed; if `Omega` is used,
  multiple cores may need more memory
  [\#77](https://github.com/r-spatial/spdep/issues/77)

- For **s2** \> 1.0.7, use indexed distances in
  [`dnearneigh()`](https://r-spatial.github.io/spdep/reference/dnearneigh.md)
  <https://github.com/r-spatial/s2/pull/162>.

## Version 1.2-1 (2022-01-05)

CRAN release: 2022-01-04

- Switching deprecated functions moved to **spatialreg** to defunct.

## Version 1.1-13 (2021-12-14)

CRAN release: 2021-12-14

- Recent changes in
  [`poly2nb()`](https://r-spatial.github.io/spdep/reference/poly2nb.md)
  had reduced and most recently (1.1-8) removed the use of `snap=` in
  finding candidate neighbours; many thanks to Matilda Brown for a clear
  and well-documented issue
  [\#65](https://github.com/r-spatial/spdep/issues/65)

- Add local Geary’s C
  [\#66](https://github.com/r-spatial/spdep/issues/66) thanks to Josiah
  Parry, discussion on further work on
  [\#68](https://github.com/r-spatial/spdep/issues/68)

- [`localmoran_perm()`](https://r-spatial.github.io/spdep/reference/localmoran.md)
  returns both look-up and folded rank p-values

## Version 1.1-12 (2021-11-09)

CRAN release: 2021-11-09

- In
  [`poly2nb()`](https://r-spatial.github.io/spdep/reference/poly2nb.md),
  reverted removal of legacy interpreted overlapping envelope code for
  sp objects that cannot be coerced to sf without **rgeos**.

- Add Fortran character handling `USE_FC_LEN_T` WRE §6.6.1.

- Checks OK with forthcoming **deldir** 1.0-0.

- Fixes [\#62](https://github.com/r-spatial/spdep/issues/62) clarifying
  [`dnearneigh()`](https://r-spatial.github.io/spdep/reference/dnearneigh.md)
  help page

## Version 1.1-11 (2021-09-07)

CRAN release: 2021-09-07

- [`knearneigh()`](https://r-spatial.github.io/spdep/reference/knearneigh.md)
  and
  [`nbdists()`](https://r-spatial.github.io/spdep/reference/nbdists.md);
  added prototype adaptation to **s2** for unprojected coordinates, used
  if [`sf_use_s2()`](https://r-spatial.github.io/sf/reference/s2.html)
  is `TRUE` which became the default for **sf** 1.0.0
  <https://github.com/r-spatial/s2/issues/125>. These are activated by
  default.

- [`dnearneigh()`](https://r-spatial.github.io/spdep/reference/dnearneigh.md)
  can choose the prototype **s2** approach if
  [`sf_use_s2()`](https://r-spatial.github.io/sf/reference/s2.html) is
  `TRUE` and `use_s2=TRUE` for unprojected coordinates; from
  <https://github.com/r-spatial/s2/issues/125> it seems that distance
  thresholds at present use brute-force rather than spatial indexing.
  Use is not activated by default.

- [`poly2nb()`](https://r-spatial.github.io/spdep/reference/poly2nb.md)
  now uses
  [`sf::st_intersects()`](https://r-spatial.github.io/sf/reference/geos_binary_pred.html)
  to find candidate neighbours unless `findInBounds=` is not NULL. With
  spatial indexing, this is very fast and scales well for large data
  sets. If
  [`sf_use_s2()`](https://r-spatial.github.io/sf/reference/s2.html) is
  `TRUE`,
  [`sf::st_intersects()`](https://r-spatial.github.io/sf/reference/geos_binary_pred.html)
  passes the geometries to
  [`s2::s2_intersects_matrix()`](https://r-spatial.github.io/s2/reference/s2_closest_feature.html),
  which also uses spatial indexing and is very fast, scaling well for
  large data sets.

- [`localmoran()`](https://r-spatial.github.io/spdep/reference/localmoran.md)
  and
  [`localmoran_perm()`](https://r-spatial.github.io/spdep/reference/localmoran.md)
  return cluster quadrants in an attribute for three splits, zeros,
  means and medians on the variable of interest and its spatial lag.

- [`localmoran_perm()`](https://r-spatial.github.io/spdep/reference/localmoran.md)
  returns the skewness and kurtosis of the permutation samples.

## Version 1.1-8 (2021-05-23)

CRAN release: 2021-05-23

- [\#55](https://github.com/r-spatial/spdep/issues/55) related to
  [\#20](https://github.com/r-spatial/spdep/issues/20) and cycling order
  in setting up grids provoked re-design of interface to
  [`cell2nb()`](https://r-spatial.github.io/spdep/reference/cell2nb.md),
  with passing of `"GridTopology"` or `"SpatialGrid"` objects as unnamed
  first or `x=` argument. Coerce `"RasterLayer"` or similar **raster**,
  **terra** or **stars** objects to **sp** class objects first if need
  be.

- In working with renewing the arguments to
  [`cell2nb()`](https://r-spatial.github.io/spdep/reference/cell2nb.md),
  it was useful to add **tinytest** support, which is now present for
  this function and may be extended to other functions for creating
  `"nb"` objects.

- [\#58](https://github.com/r-spatial/spdep/issues/58) contributed by
  Jeff Sauer and Levi Wolf (from
  <https://doi.org/10.31219/osf.io/ugkhp>) providing conditional
  standard deviates for local Moran’s I

- Error in assignment to matrix detected by CRAN check in SIDS vignette,
  section on median polish

## Version 1.1-7 (2021-04-03)

CRAN release: 2021-04-03

- Changes to continuous integration and vignettes.

- Error in `poly2nb(, queen=FALSE)` in **sf** grids (double counting of
  closed polygon start/end points),
  <https://github.com/r-spatial/spdep/issues/50>, thanks to Christopher
  Kenny.

- Adding local Moran and local G conditional permutation:
  [`localmoran_perm()`](https://r-spatial.github.io/spdep/reference/localmoran.md)
  and
  [`localG_perm()`](https://r-spatial.github.io/spdep/reference/localG.md).

- Adding
  [`nb2listwdist()`](https://r-spatial.github.io/spdep/reference/nb2listwdist.md)
  contributed by René Westerholt.

- Adding use of **sf** through GEOS to find polygon contiguity
  candidates in
  [`poly2nb()`](https://r-spatial.github.io/spdep/reference/poly2nb.md)
  if geometry count \>= 500 - uses intersections in polygon envelopes.

- [\#38](https://github.com/r-spatial/spdep/issues/38),
  [\#53](https://github.com/r-spatial/spdep/issues/53) removing
  **RANN**, adding **dbscan** suggestions for fast
  [`dnearneigh()`](https://r-spatial.github.io/spdep/reference/dnearneigh.md)
  and
  [`knearneigh()`](https://r-spatial.github.io/spdep/reference/knearneigh.md)
  via `use_kd_tree=` argument for fast planar neighbour set finding in
  2D and 3D. Affects
  [`soi.graph()`](https://r-spatial.github.io/spdep/reference/graphneigh.md)
  too, which had used **RANN**.

- [\#54](https://github.com/r-spatial/spdep/issues/54) avoid partial
  matching in `glist=` handling.

- Disambiguating **spdep** and **spatialreg** model output object class
  names prior to making **spdep** model fitting functions defunct.

## Version 1.1-5 (2020-06-29)

CRAN release: 2020-06-29

- Replacing broken geoda URLs, moving knitr to rmarkdown, work-around
  missing weights files in spData.

## Version 1.1-3 (2019-09-18)

CRAN release: 2019-09-18

- A small maintenance update to accommodate a forthcoming change in
  spData (a dataset used in an example in spdep from spData is changing
  its name; the name had involved putting “x”, “y” and “xyz” in the
  global environment through lazy loading a dataset).

## Version 1.1-2 (2019-04-05)

CRAN release: 2019-04-05

- Follow-up version of spdep with all the functions and methods
  transferred to the spatialreg package marked as deprecated, but still
  exported from spdep. Reverse dependencies passing with the released
  version still pass for me with this version.
