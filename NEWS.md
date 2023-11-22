# Version 1.3-1 (development)

* functions creating `nb` objects now warn if the object has a sub-graph count of > 1 and  `get.SubgraphOption` is `TRUE` (default `FALSE`): `complement.nb`, `diffnb`, `dnearneigh`, `droplinks`, `edit.nb`, `graph2nb`, `knn2nb`, `nb2blocknb`, `nblag`, `nblag_cumul`, `poly2nb`, `read.gal`, `read.gwt2nb`, `setdiff.nb`, `tolerance.nb`, `tri2nb`, `union.nb`

* `summary.nb`, `print.nb`, `summary.listw` and `print.listw` now report the subgraph count from `n.comp.nb` if it is more than one

* `subset.nb` now reports if the subgraph count of the neighbour object increases on subsetting

* adding a `zero.policy` attribute to functions creating `listw` objects: `nb2listw`, `sn2listw`, `mat2listw`, `nb2listwdist`. Default `zero.policy=` argument updated to use `attr(., "zero.policy")` in `summary.listw`, `print.listw`, `moran`, `moran.test`, `moran.mc`, `moran.plot`, `geary.mc`, `geary`, `geary.test`, `globalG.test`, `joincount.test`, `joincount.mc`, `joincount.multi`, `localC`, `localC_perm`, `localmoran`, `localmoran_perm`, `localG`, `localG_perm`, `lee`, `lee.test`, `lee.mc`, `lm.morantest`, `lm.LMtests`, `sp.mantel.mc`, `listw2star`, `lag.listw`, `lm.morantest`, `lm.LMtests`, `subset.listw`, `EBImoran.mc`, `LOSH`, `LOSH.mc`, `LOSH.cs`, `lm.morantest.exact`and `lm.morantest.sad`

* confusing error message in `moran.plot()` if no-neighbour cases, but `zero.policy=FALSE`

* replace `rgrass7` with `rgrass` in vignette

* fix #133 (`edit.nb` affected by not attaching `sp`)

# Version 1.2-8 (2023-02-28)

* `mat2listw()` warning if no `style=` argument given, or if `M"` is given https://github.com/r-spatial/spatialreg/issues/24, https://github.com/r-spatial/spatialreg/issues/23. 

* remaining users of `run_perm()` - `localC()`, `localmoran_bv()` and `local_joincount_uni()` get `no_repeat_in_row=` arguments.

* Address 2) in #124; `localG_perm()` and `localmoran_perm()` get `no_repeat_in_row=` arguments to use conditional permutation without replacement by sample vectors; the default implementation uses sampling with replacement, which is acceptable across simulation draws, but arguably less acceptable within draws. Feedback would be valued.

* Address 1) in #124; `localG_perm()` and `localG()` now return the same analytical standard deviates. The standard deviates from the simulated distributions are now returned in `attr(., "internals")[,"StdDev.Gi"]` from `localG_perm()`, as are p-values, etc.

* move **sp** from Depends to Imports, to reduce the visual impression that **sp** objects are required for **spdep**; **sf** objects are now preferred, but **sp** objects can be used as before, although users may need to attach **sp** expliciitly.

* fix #121 and #123; correcting returned values for `localG_perm()` when estimating the G-star measure (fix self x values and weights)

* address #120, moving documentation of `listw2U()` from `?lm.morantest` to `?nb2listw`

* addressing #119 for interpretation of `moran_bv()` results

* adding #116, René Westerholt

* fix #113, too low R version for `grDevices::hcl.colors()`

* addressing #111 by Josiah Parry

* PRs from René Westerholt, ending with #109, for the local GS measure 

# Version 1.2-7 (2022-10-01)

* #103 refactoring local joincount test by Josiah Parry

* add `hotspot` methods for `localmoran` (analytical, permutation, Saddlepoint and exact), `localC` (univariate and multivariate) and `localG` (analytical and permutation)

* #95 add `"two.sided"` to `lee.mc()`, same for `sp.mantel.mc()`, `EBImoran.mc()`, `joincount.mc()`, `geary.mc()`

* #92, #93, #94, #96, #97  contributions of prototype bivariate Moran, local  bivariate Moran and local joincount and bivariate joincount tests by Josiah Parry

* #91 `tolerance.nb()` update by F. Guillaume Blanchet

* updating coercion for **Matrix** 1.4-2

* fix ncpus issue in dontrun examples

* remove suggested packages **rgdal**, **rgeos**, **maptools**

# Version 1.2-5 (2022-08-11)

* permit use of data.frame or tibble as matrix for functions creating neighbour objects from 2D points (preferred use an object inheriting from  ``"SpatialPoints"` or `"sfc"`)

* fix #87 wrong logic in infinite weights in `nb2listwdist()`

* https://github.com/r-spatial/s2/pull/174 speeds up `dnearneigh()` for geographical coordinates without using `s2::s2_closest_edges()`.

* Adapting vignettes for absence of **rgdal** and **maptools**

# Version 1.2-4 (2022-04-18)

* added `remove.self()`, thanks to Josiah Parry #83.

* unescape underscores in help pages.

# Version 1.2-3 (2022-03-29)

* replace deprecated S-compatibility macros `DOUBLE_`

* #81 improved `dnearneigh()` help page.

* #79 remove `"htest"` class from `LOSH.mc()` output object.

* Added GA SI article to citations.

# Version 1.2-2 (2022-01-28)

* Replace `rainbow()` by `hcl.colors(..., "Set 2")` in `plot.skater()`.

* Add link to R-sig-geo thread on `EBlocal()` NaN estimates when many counts are zero on help page.

* Revise and add documentation for object returned by `localC_perm()` #68 #72 #73 #74 #75 #76

* `localmoran.sad()`, `localmoran.exact()` and `localmoran.exact.alt()` will now use multiple compute nodes if needed; if `Omega` is used, multiple cores may need more memory #77

* For **s2** > 1.0.7, use indexed distances in `dnearneigh()` https://github.com/r-spatial/s2/pull/162.

# Version 1.2-1 (2022-01-05)

* Switching deprecated functions moved to **spatialreg** to defunct.

# Version 1.1-13 (2021-12-14)

* Recent changes in `poly2nb()` had reduced and most recently (1.1-8) removed the use of `snap=` in finding candidate neighbours; many thanks to Matilda Brown for a clear and well-documented issue #65 

* Add local Geary's C #66 thanks to Josiah Parry, discussion on further work on #68

* `localmoran_perm()` returns both look-up and folded rank p-values

# Version 1.1-12 (2021-11-09)

* In `poly2nb()`, reverted removal of legacy interpreted overlapping envelope code for sp objects that cannot be coerced to sf without **rgeos**.

* Add Fortran character handling `USE_FC_LEN_T` WRE §6.6.1.

* Checks OK with forthcoming **deldir** 1.0-0.

* Fixes #62 clarifying `dnearneigh()` help page

# Version 1.1-11 (2021-09-07)

* `knearneigh()` and `nbdists()`; added prototype adaptation to **s2** for unprojected coordinates, used if `sf_use_s2()` is `TRUE` which became the default for **sf** 1.0.0 https://github.com/r-spatial/s2/issues/125. These are activated by default.

* `dnearneigh()` can choose the prototype **s2** approach if `sf_use_s2()` is `TRUE` and `use_s2=TRUE` for unprojected coordinates; from https://github.com/r-spatial/s2/issues/125 it seems that distance thresholds at present use brute-force rather than spatial indexing. Use is not activated by default.

* `poly2nb()` now uses `sf::st_intersects()` to find candidate neighbours unless `findInBounds=` is not NULL. With spatial indexing, this is very fast and scales well for large data sets. If `sf_use_s2()` is `TRUE`, `sf::st_intersects()` passes the geometries to `s2::s2_intersects_matrix()`, which also uses spatial indexing and is very fast, scaling well for large data sets.

* `localmoran()` and `localmoran_perm()` return cluster quadrants in an attribute for three splits, zeros, means and medians on the variable of interest and its spatial lag.

* `localmoran_perm()` returns the skewness and kurtosis of the permutation samples.


# Version 1.1-8 (2021-05-23)

* #55 related to #20 and cycling order in setting up grids provoked re-design of interface to `cell2nb()`, with passing of `"GridTopology"` or `"SpatialGrid"` objects as unnamed first or `x=` argument. Coerce `"RasterLayer"` or similar **raster**, **terra** or **stars** objects to **sp** class objects first if need be.

* In working with renewing the arguments to `cell2nb()`, it was useful to add **tinytest** support, which is now present for this function and may be extended to other functions for creating `"nb"` objects.

* #58 contributed by Jeff Sauer and Levi Wolf (from https://doi.org/10.31219/osf.io/ugkhp) providing conditional standard deviates for local Moran's I

* Error in assignment to matrix detected by CRAN check in SIDS vignette, section on median polish

# Version 1.1-7 (2021-04-03)

* Changes to continuous integration and vignettes.

* Error in `poly2nb(, queen=FALSE)` in **sf** grids (double counting of closed polygon start/end points), https://github.com/r-spatial/spdep/issues/50, thanks to Christopher Kenny.

* Adding local Moran and local G conditional permutation: `localmoran_perm()` and `localG_perm()`.

* Adding `nb2listwdist()` contributed by René Westerholt.

* Adding use of **sf** through GEOS to find polygon contiguity candidates in `poly2nb()` if geometry count >= 500 - uses intersections in polygon envelopes.

* #38, #53 removing **RANN**, adding **dbscan** suggestions for fast `dnearneigh()` and `knearneigh()` via `use_kd_tree=` argument for fast planar neighbour set finding in 2D and 3D. Affects `soi.graph()` too, which had used **RANN**.

* #54 avoid partial matching in `glist=` handling.

* Disambiguating **spdep** and **spatialreg** model output object class names prior to making **spdep** model fitting functions defunct.

# Version 1.1-5 (2020-06-29)

* Replacing broken geoda URLs, moving knitr to rmarkdown, work-around missing weights files in spData.


# Version 1.1-3 (2019-09-18)

* A small maintenance update to accommodate a forthcoming change in spData (a dataset used in an example in spdep from spData is changing its name;   the name had involved putting "x", "y" and "xyz" in the global environment through lazy loading a dataset).


# Version 1.1-2 (2019-04-05)

* Follow-up version of spdep with all the functions and 
  methods transferred to the spatialreg package marked 
  as deprecated, but still exported from spdep. Reverse 
  dependencies passing with the released version still 
  pass for me with this version.
