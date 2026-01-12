# Read and write spatial neighbour files

The "gwt" functions read and write GeoDa GWT files (the example file
baltk4.GWT was downloaded from the site given in the reference), and the
"dat" functions read and write Matlab sparse matrix files as used by
James LeSage's Spatial Econometrics Toolbox (the example file wmat.dat
was downloaded from the site given in the reference). The body of the
files after any headers should have three columns separated by white
space, and the third column must be numeric in the locale of the reading
platform (correct decimal separator).

## Usage

``` r
read.gwt2nb(file, region.id=NULL)
write.sn2gwt(sn, file, shpfile=NULL, ind=NULL, useInd=FALSE, legacy=FALSE)
read.dat2listw(file)
write.sn2dat(sn, file)
read.swmdbf2listw(fn, region.id=NULL, style=NULL, zero.policy=NULL)
read_swm_dbf(fn)
write.swmdbf(listw, file, ind, region.id = attr(listw, "region.id"))
write_swm_dbf(listw, file, ind, region.id = attr(listw, "region.id"))
write.sn2DBF(sn, file)
```

## Arguments

- file, fn:

  name of file with weights data

- region.id:

  a character vector of region IDs - for ArcGIS SWM DBFs, the values
  must be character integers (only numbers not starting with zero)

- listw:

  a `listw` object

- sn:

  a `spatial.neighbour` object

- shpfile:

  character string: if not given Shapefile name taken from GWT file for
  this dataset

- ind:

  character string: region id indicator field name

- useInd:

  default FALSE, if TRUE, write `region.id` attribute ID key tags to
  output file (use in OpenGeoDa will depend on the shapefile having the
  field named in the `ind` argument matching the exported tags)

- legacy:

  default FALSE; if TRUE, header has single field with number of
  observations only

- style:

  default NULL, missing, set to "M" and warning given; if not "M",
  passed to
  [`nb2listw`](https://r-spatial.github.io/spdep/reference/nb2listw.md)
  to re-build the object

- zero.policy:

  default NULL, use global option value; if FALSE stop with error for
  any empty neighbour sets, if TRUE permit the weights list to be formed
  with zero-length weights vectors

## Details

Attempts to honour the region.id argument given when reading GWT and
SWM/DBF files. If the region IDs given in `region.id=` do not match the
origins or destinations in the GWT file, an error or warning will be
thrown, which should be considered carefully. `read_swm_dbf` is a
simplified interface to `read.swmdbf2listw` When no-neighbour
observations are present, it is advised that `region.id=` be given in
`read.swmdbf2listw`; while the function should read correctly if the
minimum and maximum IDs are present as observations with neighbours
without `region.id=` set, reading using `read.swmdbf2listw` will fail
when the minimum or maximum ID observations have no neighbours and
`region.id=` is not given.

## Value

`read.gwt2nb` returns a neighbour "nb" object with the generalised
weights stored as a list element called "dlist" of the "GeoDa"
attribute; `read.swmdbf2listw` returns a "listw" object read from a DBF
file exported from an ArcGIS SWM object.

## References

Luc Anselin (2003) *GeoDa 0.9 User's Guide*, pp. 80–81, Spatial Analysis
Laboratory, Department of Agricultural and Consumer Economics,
University of Illinois, Urbana-Champaign,
<http://geodacenter.github.io/docs/geoda093.pdf>; also material formerly
at spatial-econometrics.com/data/contents.html

## Author

Roger Bivand <Roger.Bivand@nhh.no>

## See also

[`read.gal`](https://r-spatial.github.io/spdep/reference/read.gal.md)

## Examples

``` r
data(baltimore, package="spData")
STATION <- baltimore$STATION
gwt1 <- read.gwt2nb(system.file("weights/baltk4.GWT", package="spData")[1],
 STATION)
#> Warning: 102, 115, 208 are not destinations
cat(paste("Neighbours list symmetry;", is.symmetric.nb(gwt1, FALSE, TRUE),
 "\n"))
#> Neighbours list symmetry; FALSE 
listw1 <- nb2listw(gwt1, style="B", glist=attr(gwt1, "GeoDa")$dist)
tmpGWT <- tempfile()
write.sn2gwt(listw2sn(listw1), tmpGWT)
gwt2 <- read.gwt2nb(tmpGWT, STATION)
#> Warning: 102, 115, 208 are not destinations
cat(paste("Neighbours list symmetry;", is.symmetric.nb(gwt2, FALSE, TRUE),
 "\n"))
#> Neighbours list symmetry; FALSE 
diffnb(gwt1, gwt2)
#> Warning: neighbour object has 211 sub-graphs
#> Neighbour list object:
#> Number of regions: 211 
#> Number of nonzero links: 0 
#> Percentage nonzero weights: 0 
#> Average number of links: 0 
#> 211 regions with no links:
#> 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
#> 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38,
#> 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56,
#> 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74,
#> 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92,
#> 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107,
#> 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121,
#> 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135,
#> 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149,
#> 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163,
#> 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177,
#> 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
#> 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205,
#> 206, 207, 208, 209, 210, 211
#> 211 disjoint connected subgraphs
data(oldcol)
tmpMAT <- tempfile()
COL.W <- nb2listw(COL.nb)
write.sn2dat(listw2sn(COL.W), tmpMAT)
listwmat1 <- read.dat2listw(tmpMAT)
#> Warning: style is M (missing); style should be set to a valid value
diffnb(listwmat1$neighbours, COL.nb, verbose=TRUE)
#> Warning: region.id differ; using ids of first list
#> Warning: neighbour object has 49 sub-graphs
#> Neighbour list object:
#> Number of regions: 49 
#> Number of nonzero links: 0 
#> Percentage nonzero weights: 0 
#> Average number of links: 0 
#> 49 regions with no links:
#> 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
#> 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38,
#> 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49
#> 49 disjoint connected subgraphs
listwmat2 <- read.dat2listw(system.file("etc/weights/wmat.dat", 
 package="spdep")[1])
#> Warning: style is M (missing); style should be set to a valid value
diffnb(listwmat1$neighbours, listwmat2$neighbours, verbose=TRUE)
#> Warning: neighbour object has 49 sub-graphs
#> Neighbour list object:
#> Number of regions: 49 
#> Number of nonzero links: 0 
#> Percentage nonzero weights: 0 
#> Average number of links: 0 
#> 49 regions with no links:
#> 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
#> 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38,
#> 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49
#> 49 disjoint connected subgraphs
run <- FALSE
if (require("foreign", quietly=TRUE)) run <- TRUE
if (run) {
nc_sf <- sf::st_read(system.file("gpkg/nc.gpkg", package="sf")[1])
nc_sf$UniqueID <- 1:nrow(nc_sf)
fn <- system.file("etc/misc/nc_contiguity_unique_id.dbf", package="spdep")[1]
nc1 <- read.swmdbf2listw(fn, style="B")
nc1
}
#> Reading layer `nc.gpkg' from data source `/home/rsb/lib/r_libs/sf/gpkg/nc.gpkg' using driver `GPKG'
#> Simple feature collection with 100 features and 14 fields
#> Geometry type: MULTIPOLYGON
#> Dimension:     XY
#> Bounding box:  xmin: -84.32385 ymin: 33.88199 xmax: -75.45698 ymax: 36.58965
#> Geodetic CRS:  NAD27
#> Warning: region.id not given, c(MYID, NID) range is 1:100
#> Characteristics of weights list object:
#> Neighbour list object:
#> Number of regions: 100 
#> Number of nonzero links: 490 
#> Percentage nonzero weights: 4.9 
#> Average number of links: 4.9 
#> 
#> Weights style: B 
#> Weights constants summary:
#>     n    nn  S0  S1    S2
#> B 100 10000 490 980 10696
if (run) {
nc1a <- read.swmdbf2listw(fn, region.id=as.character(nc_sf$UniqueID),
 style="B")
all.equal(nc1, nc1a)
}
#> [1] TRUE
if (run) {
fn <- system.file("etc/misc/nc_contiguity_unique_id_islands.dbf",
 package="spdep")[1]
try(nc1i <- read.swmdbf2listw(fn, style="B"))
nc1i <- read.swmdbf2listw(fn, style="B", zero.policy=TRUE)
nc1ia <- read.swmdbf2listw(fn, region.id=as.character(nc_sf$UniqueID),
 style="B", zero.policy=TRUE)
nc1ia
}
#> Warning: region.id not given, c(MYID, NID) range is 1:100
#> Warning: 81 is not an origin
#> Error in read.swmdbf2listw(fn, style = "B") : 
#>   Error in sn2listw(df1, style = style, zero.policy = zero.policy) : 
#>   no-neighbour observations found, set zero.policy to TRUE
#> 
#> Warning: region.id not given, c(MYID, NID) range is 1:100
#> Warning: no-neighbour observations found, set zero.policy to TRUE;
#> this warning will soon become an error
#> Warning: no-neighbour observations found, set zero.policy to TRUE;
#> this warning will soon become an error
#> Characteristics of weights list object:
#> Neighbour list object:
#> Number of regions: 100 
#> Number of nonzero links: 487 
#> Percentage nonzero weights: 4.87 
#> Average number of links: 4.87 
#> 1 region with no links:
#> 81
#> Non-symmetric neighbours list
#> 
#> Weights style: B 
#> Weights constants summary:
#>    n   nn  S0  S1    S2
#> B 99 9801 487 971 10632
if (run) {
all.equal(nc1i, nc1ia)
}
#> [1] TRUE
if (run) {
GDAL37 <- numeric_version(unname(sf::sf_extSoftVersion()["GDAL"]), strict=FALSE)
(GDAL37 <- ifelse(is.na(GDAL37), FALSE, GDAL37 >= "3.7.0"))
file <- "etc/shapes/california.gpkg.zip"
zipfile <- system.file(file, package="spdep")
if (GDAL37) {
    cal <- st_read(zipfile)
} else {
    td <- tempdir()
    bn <- sub(".zip", "", basename(file), fixed=TRUE)
    target <- unzip(zipfile, files=bn, exdir=td)
    cal <- st_read(target)
}
fn <- system.file("etc/misc/contiguity_myid.dbf", package="spdep")[1]
cal1 <- read.swmdbf2listw(fn, style="B")
cal1a <- read.swmdbf2listw(fn, region.id=as.character(cal$MYID), style="B")
all.equal(cal1, cal1a)
}
#> Reading layer `california' from data source 
#>   `/tmp/Rtmp6xgijc/temp_libpath8c3932b2ce99/spdep/etc/shapes/california.gpkg.zip' 
#>   using driver `GPKG'
#> Simple feature collection with 58 features and 2 fields
#> Geometry type: MULTIPOLYGON
#> Dimension:     XY
#> Bounding box:  xmin: -13847330 ymin: 3810869 xmax: -12704360 ymax: 5132708
#> Projected CRS: World_Mercator
#> Warning: region.id not given, c(MYID, NID) range is 158:215
#> [1] TRUE
if (run) {
fn <- system.file("etc/misc/contiguity_unique_id.dbf", package="spdep")[1]
cal2 <- read.swmdbf2listw(fn, style="B")
cal2a <- read.swmdbf2listw(fn, region.id=as.character(cal$UniqueID), style="B")
all.equal(cal2, cal2a)
}
#> Warning: region.id not given, c(MYID, NID) range is 1:58
#> [1] TRUE
if (run) {
fn <- system.file("etc/misc/contiguity_unique_id_islands.dbf", package="spdep")[1]
try(cal3i <- read.swmdbf2listw(fn, style="B"))
cal3i <- read.swmdbf2listw(fn, style="B", zero.policy=TRUE)
cal3ia <- read.swmdbf2listw(fn, region.id=as.character(cal$UniqueID), style="B", zero.policy=TRUE)
all.equal(cal3i, cal3ia)
}
#> Warning: region.id not given, c(MYID, NID) range is 1:58
#> Warning: 21, 38 are not origins
#> Error in read.swmdbf2listw(fn, style = "B") : 
#>   Error in sn2listw(df1, style = style, zero.policy = zero.policy) : 
#>   no-neighbour observations found, set zero.policy to TRUE
#> 
#> Warning: region.id not given, c(MYID, NID) range is 1:58
#> Warning: no-neighbour observations found, set zero.policy to TRUE;
#> this warning will soon become an error
#> Warning: no-neighbour observations found, set zero.policy to TRUE;
#> this warning will soon become an error
#> [1] TRUE
if (run) {
cal1a_1n_nb <- cal1a$neighbours
cal1a_1n_nb <- droplinks(cal1a_1n_nb, drop=c("158", "180", "215"), sym=TRUE)
cal1a_1n <- nb2listw(cal1a_1n_nb, style="B", zero.policy=TRUE)
cal1a_1n_sn <- listw2sn(cal1a_1n)
file <- tempfile(fileext=".dbf")
write.sn2DBF(cal1a_1n_sn, file)
cal1a_1n_rt <- read.swmdbf2listw(file, region.id=as.character(cal$MYID),
 style="B", zero.policy=TRUE)
all.equal(cal1a_1n$neighbours, cal1a_1n_rt$neighbours)
}
#> Warning: some observations have no neighbours
#> Warning: neighbour object has 4 sub-graphs
#> Warning: no-neighbour observations found, set zero.policy to TRUE;
#> this warning will soon become an error
#> Warning: neighbour object has 4 sub-graphs
#> [1] TRUE
if (run) {
all.equal(cal1a_1n$weights, cal1a_1n_rt$weights, check.attributes=FALSE)
}
#> [1] TRUE
if (run) {
cal1_1n_rt <- read.swmdbf2listw(file, style="B", zero.policy=TRUE)
all(isTRUE(all.equal(cal1a_1n$neighbours, cal1_1n_rt$neighbours)))
}
#> Warning: region.id not given, c(MYID, NID) range is 159:214
#> Warning: no-neighbour observations found, set zero.policy to TRUE;
#> this warning will soon become an error
#> Warning: neighbour object has 2 sub-graphs
#> [1] FALSE
if (run) {
all(isTRUE(all.equal(cal1a_1n$weights, cal1_1n_rt$weights)))
}
#> [1] FALSE
```
