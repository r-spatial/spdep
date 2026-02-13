# Creating Neighbours using sf objects

## Creating Neighbours using sf objects

### Introduction

This vignette tracks the legacy nb vignette, which was based on part of
the first (2008) edition of ASDAR. It adds hints to the code in the nb
vignette to use the sf vector representation instead of the sp vector
representation to create neighbour objects. .

### Summary

This is a summary of the results below:

- In general, if you need to reproduce results from using `"Spatial"`
  objects in **spdep**, coerce sf objects to sp objects before
  constructing neighbour objects (particularly if polygon centroids are
  used for point representation).

- However, for new work, you should use `"sf"` objects read in using
  **sf**.

- From **spdep** 1.1-7, a number of steps have been taken to choose more
  efficient approaches especially for larger data sets, using functions
  in **sf** and the approximate nearest neighbour (ANN) implementation
  in **dbscan** rather than **RANN**.

### Data set

We’ll use the whole NY 8 county set of boundaries, as they challenge the
implementations more than just the Syracuse subset. The description of
the input geometries from ADSAR is: New York leukemia: used and
documented extensively in Waller and Gotway (2004) and with data
formerly made available in Chap. 9 of
`http://web1.sph.emory.edu/users/lwaller/ch9index.htm`; the data import
process is described in the help file of NY_data in **spData**;
geometries downloaded from the CIESIN server formerly at
sedac.ciesin.columbia.edu/data/set/acrp-boundary-1992/data-download, now
possibly
earthdata.nasa.gov/data/catalog/sedac-ciesin-sedac-acrp-1992-bf-1.00,
file /pub/census/usa/tiger/ny/bna_st/t8_36.zip, and extensively edited;
a zip archive NY_data.zip of shapefiles and a GAL format neighbours list
is on the book website. Further, the zipfile is now at a new location
requiring login. The object listw_NY is directly imported from
nyadjwts.dbf on the Waller & Gotway (2004) chapter 9 website.

The version of the New York 8 counties geometries used in ASDAR and
included as a shapefile in spdep was converted from the original BNA
file using an external utility program to convert to MapInfo format and
converted on from there using GDAL 1.4.1 (the OGR BNA driver was not
then available; it entered OGR at 1.5.0, release at the end of 2007),
and contains invalid geometries. What was found in mid-2007 was that
included villages were in/excluded by in-out umbilical cords to the
boundary of the enclosing tract, when the underlying BNA file was first
converted to MapInfo (holes could not exist then).

Here we will use a GPKG file created as follows (rgdal could also be
used with the same output; GDAL prior to GDAL 3.3 built with GEOS, so
the BNA vector driver will use geometry tests: The BNA driver supports
reading of polygons with holes or lakes. It determines what is a hole or
a lake only from geometrical analysis (inclusion, non-intersection
tests) and ignores completely the notion of polygon winding (whether the
polygon edges are described clockwise or counter-clockwise). GDAL must
be built with GEOS enabled to make geometry test work.):

``` r
library(sf)
sf_bna <- st_read("t8_36.bna", stringsAsFactors=FALSE)
table(st_is_valid(sf_bna))
sf_bna$AREAKEY <- gsub("\\.", "", sf_bna$Primary.ID)
data(NY_data, package="spData")
key <- as.character(nydata$AREAKEY)
sf_bna1 <- sf_bna[match(key, sf_bna$AREAKEY), c("AREAKEY")]
sf_bna2 <- merge(sf_bna1, nydata, by="AREAKEY")
sf_bna2_utm18 <- st_transform(sf_bna2, "+proj=utm +zone=18 +datum=NAD27")
st_write(sf_bna2_utm18, "NY8_bna_utm18.gpkg")
```

### nb and listw objects (copied from the nb_igraph vignette)

Since the **spdep** package was created, *spatial weights* objects have
been constructed as lists with three components and a few attributes, in
old-style class `listw` objects. The first component of a `listw` object
is an `nb` object, a list of `n` integer vectors, with at least a
character vector `region.id` attribute with `n` unique values (like the
`row.names` of a `data.frame` object); `n` is the number of spatial
entities. Component `i` of this list contains the integer identifiers of
the neighbours of `i` as a sorted vector with no duplication and values
in `1:n`; if `i` has no neighbours, the component is a vector of length
`1` with value `0L`. The `nb` object may contain an attribute indicating
whether it is symmetric or not, that is whether `i` is a neighbour of
`j` implies that `j` is a neighbour of `i`. Some neighbour definitions
are symmetric by construction, such as contiguities or distance
thresholds, others are asymmetric, such as `k`-nearest neighbours. The
`nb` object redundantly stores both `i`-`j` and `j`-`i` links.

The second component of a `listw` object is a list of `n` numeric
vectors, each of the same length as the corresponding non-zero vectors
in the `nb`object. These give the values of the spatial weights for each
`i`-`j` neighbour pair. It is often the case that while the neighbours
are symmetric by construction, the weights are not, as for example when
weights are *row-standardised* by dividing each row of input weights by
the count of neighbours or cardinality of the neighbour set of `i`. In
the `nb2listw`function, it is also possible to pass through general
weights, such as inverse distances, shares of boundary lengths and so
on.

The third component of a `listw` object records the `style` of the
weights as a character code, with `"B"` for binary weights taking values
zero or one (only one is recorded), `"W"` for row-standardised weights,
and so on. In order to subset `listw` objects, knowledge of the `style`
may be necessary.

### Comparison of sp and sf approaches

First some housekeeping and setup to permit this vignette to be built
when packages are missing or out-of-date:

``` r
if (!suppressPackageStartupMessages(require(sf, quietly=TRUE))) {
  message("install the sf package")
  dothis <- FALSE
}
if (dothis) sf_extSoftVersion()
```

    ##           GEOS           GDAL         proj.4 GDAL_with_GEOS     USE_PROJ_H 
    ##       "3.14.1"       "3.12.2"        "9.7.1"         "true"         "true" 
    ##           PROJ 
    ##        "9.7.1"

Let us read the GPKG file with valid geometries in to ‘sf’ and ‘sp’
objects:

``` r
NY8_sf <- st_read(system.file("shapes/NY8_bna_utm18.gpkg", package="spData"), quiet=TRUE)
table(st_is_valid(NY8_sf))
```

    ## 
    ## TRUE 
    ##  281

#### Contiguity neighbours for polygon support

Here we first generate a queen contiguity nb object using the legacy
spdep approach. This first either uses a pre-computed list of vectors of
probable neighbours or finds intersecting bounding boxes internally.
Then the points on the boundaries of each set of polygons making up an
observation are checked for a distance less than snap to any of the
points of the set of polygons making up an observation included in the
set of candidate neighbours. Because contiguity is symmetric, only i to
j contiguities are tested. A queen contiguity is found as soon as one
point matches, a rook contiguity as soon as two points match:

``` r
suppressPackageStartupMessages(library(spdep))
reps <- 10
eps <- sqrt(.Machine$double.eps)
system.time(for(i in 1:reps) NY8_sf_1_nb <- poly2nb(NY8_sf, queen=TRUE, snap=eps))/reps
```

    ##    user  system elapsed 
    ##  0.1895  0.0100  0.2003

Using spatial indices to check intersection of polygons is much faster
than the legacy method in poly2nb. From **spdep** 1.1-7, use is made of
GEOS through **sf** to find candidate neighbours when `foundInBox=NULL`,
the default value. Because contiguity is symmetric by definition,
`foundInBox=` only requires intersections for higher indices, leading to
a slight overhead to remove duplicates, as
[`st_intersects()`](https://r-spatial.github.io/sf/reference/geos_binary_pred.html)
reports both `i j` ans `j i` relationships. As
[`st_intersects()`](https://r-spatial.github.io/sf/reference/geos_binary_pred.html)
does not report whether neighbours are `queen` or `rook`, a further step
is needed to distinguish the two cases.

``` r
NY8_sf_1_nb
```

    ## Neighbour list object:
    ## Number of regions: 281 
    ## Number of nonzero links: 1632 
    ## Percentage nonzero weights: 2.066843 
    ## Average number of links: 5.807829

spdep::poly2nb uses two heuristics, first to find candidate neighbours
from intersecting polygons
([`st_intersects()`](https://r-spatial.github.io/sf/reference/geos_binary_pred.html)),
and second to use the symmetry of the relationship to halve the number
of remaining tests. This means that performance is linear in n, but with
overhead for identifying candidates, and back-filling symmetric
neighbours. Further,
[`spdep::poly2nb()`](https://r-spatial.github.io/spdep/reference/poly2nb.md)
stops searching for queen contiguity as soon as the first neighbour
point is found within snap distance (if not identical, which is tested
first); a second neighbour point indicates rook contiguities. For
details of alternatives for spherical geometries, see section
@ref(spher-poly2nb) below.

#### Contiguity neighbours from invalid polygons

Next, we explore a further possible source of differences in neighbour
object reproduction, using the original version of the tract boundaries
used in ASDAR, but with some invalid geometries as mentioned earlier
(`NY8_utm18.gpkg` was created from the original ESRI Shapefile used in
ASDAR):

``` r
if (packageVersion("spData") >= "2.3.2") {
    NY8_sf_old <- sf::st_read(system.file("shapes/NY8_utm18.gpkg", package="spData"))
} else {
    NY8_sf_old <- sf::st_read(system.file("shapes/NY8_utm18.shp", package="spData"))
}
```

    ## Reading layer `NY8_utm18' from data source 
    ##   `/home/rsb/lib/r_libs/spData/shapes/NY8_utm18.gpkg' using driver `GPKG'
    ## Simple feature collection with 281 features and 17 fields
    ## Geometry type: POLYGON
    ## Dimension:     XY
    ## Bounding box:  xmin: 358241.9 ymin: 4649755 xmax: 480393.1 ymax: 4808545
    ## Projected CRS: WGS 84 / UTM zone 18N

``` r
table(st_is_valid(NY8_sf_old))
```

    ## 
    ## FALSE  TRUE 
    ##     5   276

We can see that there are a number of differences between the neighbour
sets derived from the fully valid geometries and the older partly
invalid set:

``` r
try(NY8_sf_old_1_nb <- poly2nb(NY8_sf_old), silent = TRUE)
all.equal(NY8_sf_old_1_nb, NY8_sf_1_nb, check.attributes=FALSE)
```

    ## [1] "Component 57: Numeric: lengths (4, 5) differ" 
    ## [2] "Component 58: Numeric: lengths (5, 6) differ" 
    ## [3] "Component 66: Numeric: lengths (7, 11) differ"
    ## [4] "Component 73: Numeric: lengths (4, 5) differ" 
    ## [5] "Component 260: Numeric: lengths (8, 9) differ"

Using
[`st_make_valid()`](https://r-spatial.github.io/sf/reference/valid.html)
to make the geometries valid:

``` r
NY8_sf_old_val <- st_make_valid(NY8_sf_old, dist=0)
table(st_is_valid(NY8_sf_old_val))
```

    ## 
    ## TRUE 
    ##  281

we also see that the geometry type of the geometry column changes:

``` r
class(st_geometry(NY8_sf_old))
```

    ## [1] "sfc_POLYGON" "sfc"

``` r
class(st_geometry(NY8_sf_old_val))
```

    ## [1] "sfc_GEOMETRY" "sfc"

and checking the `"sfg"` objects, two now have objects of different
topological dimensions.

``` r
table(sapply(st_geometry(NY8_sf_old_val), function(x) class(x)[[2]]))
```

    ## 
    ## MULTIPOLYGON      POLYGON 
    ##            3          278

This can be remedied using
[`st_collection_extract()`](https://r-spatial.github.io/sf/reference/st_collection_extract.html)
to get the polygon objects:

``` r
NY8_sf_old_val <- st_collection_extract(NY8_sf_old_val, "POLYGON")
table(sapply(st_geometry(NY8_sf_old_val), function(x) class(x)[[2]]))
```

    ## 
    ## MULTIPOLYGON 
    ##          281

However, in making the geometries valid, we change the geometries, so
the new sets of neighbours still differ from those made with the valid
geometries in the same ways as before imposing validity:

``` r
try(NY8_sf_old_1_nb_val <- poly2nb(NY8_sf_old_val), silent = TRUE)
all.equal(NY8_sf_old_1_nb_val, NY8_sf_1_nb, check.attributes=FALSE)
```

    ## [1] "Component 57: Numeric: lengths (4, 5) differ" 
    ## [2] "Component 58: Numeric: lengths (5, 6) differ" 
    ## [3] "Component 66: Numeric: lengths (7, 11) differ"
    ## [4] "Component 73: Numeric: lengths (4, 5) differ" 
    ## [5] "Component 260: Numeric: lengths (8, 9) differ"

The neighbour sets are the same for the old boundaries with or without
imposing validity:

``` r
all.equal(NY8_sf_old_1_nb_val, NY8_sf_old_1_nb, check.attributes=FALSE)
```

    ## [1] TRUE

### Planar point-based neighbours

#### Finding points for polygon objects

[`knearneigh()`](https://r-spatial.github.io/spdep/reference/knearneigh.md)
and
[`dnearneigh()`](https://r-spatial.github.io/spdep/reference/dnearneigh.md)
require point support, so decisions must be taken about how to place the
point in the areal object. We can use
[`st_centroid()`](https://r-spatial.github.io/sf/reference/geos_unary.html)
to get the centroids, using the `of_largest_polygon=TRUE` argument to
make sure that the centroid is that of the largest polygon id the
observation is made up of more than one external ring:

``` r
NY8_ct_sf <- st_centroid(st_geometry(NY8_sf), of_largest_polygon=TRUE)
```

or
[`st_point_on_surface()`](https://r-spatial.github.io/sf/reference/geos_unary.html)
which guarantees that the point will fall on the surface of a member
polygon:

``` r
NY8_pos_sf <- st_point_on_surface(st_geometry(NY8_sf))
```

or indeed taking the centre of the largest inscribed circle (the
function returns a radius line segment, so we choose the central point,
not the point on the circle):

``` r
if (unname(sf_extSoftVersion()["GEOS"] >= "3.9.0")) 
    NY8_cic_sf <- st_cast(st_inscribed_circle(st_geometry(NY8_sf), nQuadSegs=0), "POINT")[(1:(2*nrow(NY8_sf)) %% 2) != 0]
```

We need to check whether coordinates are planar or not:

``` r
st_is_longlat(NY8_ct_sf)
```

    ## [1] FALSE

#### Graph-based neighbours

From this, we can check the graph-based neighbours (planar coordinates
only):

``` r
suppressPackageStartupMessages(require(deldir))
NY84_nb <- tri2nb(NY8_ct_sf)
if (require(dbscan, quietly=TRUE)) {
  NY85_nb <- graph2nb(soi.graph(NY84_nb, NY8_ct_sf))
} else NY85_nb <- NULL
```

    ## 
    ## Attaching package: 'dbscan'

    ## The following object is masked from 'package:stats':
    ## 
    ##     as.dendrogram

    ## Warning in graph2nb(soi.graph(NY84_nb, NY8_ct_sf)): neighbour object has 2
    ## sub-graphs

``` r
NY86_nb <- graph2nb(gabrielneigh(NY8_ct_sf))
NY87_nb <- graph2nb(relativeneigh(NY8_ct_sf))
```

#### K-nearest neighbours

K-nearest neighbours use the coordinate matrices, and can handle Great
Circle distances, but this is not demonstrated here, as the data set
used is planar, in which case
[`dbscan::kNN()`](https://rdrr.io/pkg/dbscan/man/kNN.html) in 2D or 3D
building a kd-tree is used:

``` r
system.time(for (i in 1:reps) suppressWarnings(NY88_nb_sf <- knn2nb(knearneigh(NY8_ct_sf, k=1))))/reps
```

    ##    user  system elapsed 
    ##  0.0252  0.0012  0.0266

Legacy code may be used omitting the kd-tree:

``` r
system.time(for (i in 1:reps) suppressWarnings(NY89_nb_sf <- knn2nb(knearneigh(NY8_ct_sf, k=1, use_kd_tree=FALSE))))/reps
```

    ##    user  system elapsed 
    ##  0.0248  0.0015  0.0264

#### Distance neighbours

Distance neighbours need a threshold - `nbdists` shows the maximum
distance to first nearest neighbour:

``` r
dsts <- unlist(nbdists(NY88_nb_sf, NY8_ct_sf))
summary(dsts)
```

    ##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    ##    82.85   912.85  1801.11  3441.04  4461.26 17033.11

``` r
max_1nn <- max(dsts)
```

`dnearneigh` can also handle Great Circle distances, but this is not
demonstrated here, as the data set used is planar:

``` r
system.time(for (i in 1:reps) suppressWarnings(NY810_nb <- dnearneigh(NY8_ct_sf, d1=0, d2=0.75*max_1nn)))/reps
```

    ##    user  system elapsed 
    ##  0.0417  0.0015  0.0433

By default, the function uses
[`dbscan::frNN()`](https://rdrr.io/pkg/dbscan/man/frNN.html) to build a
kd-tree in 2D or 3D which is then used to find distance neighbours. For
small n, the argument `use_kd_tree=FALSE` may speed up computation a
little by reverting to legacy code not building a kd-tree first, but in
general the differences are so small that the user will not notice:

``` r
system.time(for (i in 1:reps) suppressWarnings(NY811_nb <- dnearneigh(NY8_ct_sf, d1=0, d2=0.75*max_1nn, use_kd_tree=FALSE)))/reps
```

    ##    user  system elapsed 
    ##  0.0251  0.0015  0.0267

### Spherical point-based neighbours

Spherical point-based neighbours may be found using Great Circle
distances. These have been used for many years, but the switch of **sf**
1.0-0 to use **s2** by default has opened up new opportunities where
spatial indexing on the sphere may help.

``` r
pts_ll <- st_transform(NY8_ct_sf, "OGC:CRS84")
st_is_longlat(pts_ll)
```

    ## [1] TRUE

#### K-nearest neighbours

If the input geometries are in geographical coordinates, and
[`sf_use_s2()`](https://r-spatial.github.io/sf/reference/s2.html) is
`TRUE`,
[`knearneigh()`](https://r-spatial.github.io/spdep/reference/knearneigh.md)
will use spatially indexed points and
[`s2::s2_closest_edges()`](https://r-spatial.github.io/s2/reference/s2_closest_feature.html)
(see
<https://github.com/r-spatial/s2/issues/125#issuecomment-860107442>)

``` r
(old_use_s2 <- sf_use_s2())
```

    ## [1] TRUE

and performs well with also with larger data sets:

``` r
sf_use_s2(TRUE)
system.time(for (i in 1:reps) pts_ll1_nb <- knn2nb(knearneigh(pts_ll, k=6)))/reps
```

    ##    user  system elapsed 
    ##  0.0344  0.0002  0.0347

For this smaller data set, the legacy approach without spatial indexing
is adequate, but slows down as the number of observations increases:

``` r
sf_use_s2(FALSE)
```

    ## Spherical geometry (s2) switched off

``` r
system.time(for (i in 1:reps) pts_ll2_nb <- knn2nb(knearneigh(pts_ll, k=6)))/reps
```

    ##    user  system elapsed 
    ##   0.022   0.000   0.022

The WGS84 ellipsoid Great Circle distances differ a very little from the
**s2** spherical distances, yielding output that here diverges for two
tract centroids:

``` r
all.equal(pts_ll1_nb, pts_ll2_nb, check.attributes=FALSE)
```

    ## [1] "Component 52: Mean relative difference: 1.466667"  
    ## [2] "Component 124: Mean relative difference: 0.0251046"

``` r
pts_ll1_nb[[52]]
```

    ## [1] 15 38 48 49 50 53

``` r
pts_ll2_nb[[52]]
```

    ## [1] 37 38 48 49 50 53

``` r
pts_ll1_nb[[124]]
```

    ## [1] 117 122 123 125 133 134

``` r
pts_ll2_nb[[124]]
```

    ## [1] 116 117 123 125 133 134

``` r
sf_use_s2(old_use_s2)
```

    ## Spherical geometry (s2) switched on

#### Distance neighbours

Distance neighbours are more problematic. While
[`nbdists()`](https://r-spatial.github.io/spdep/reference/nbdists.md)
works well with **s2** spherical coordinates, none of the tried
adaptations for
[`dnearneigh()`](https://r-spatial.github.io/spdep/reference/dnearneigh.md)
work adequately yet. An argument `use_s2=` is set to TRUE if **s2** \>
1.0-7, using
[`s2::s2_closest_edges()`](https://r-spatial.github.io/s2/reference/s2_closest_feature.html)
or the legacy brute-force approach, then only calculating distances from
`i` to `j` and copying those to `j` to `i` through symmetry. The
distance metric is alway `"km"`.

``` r
max_1nn_ll <- max(unlist(nbdists(knn2nb(knearneigh(pts_ll, k=1)), pts_ll)))
```

    ## Warning in knn2nb(knearneigh(pts_ll, k = 1)): neighbour object has 62
    ## sub-graphs

``` r
args(dnearneigh)
```

    ## function (x, d1, d2, row.names = NULL, longlat = NULL, bounds = c("GE", 
    ##     "LE"), use_kd_tree = TRUE, symtest = FALSE, use_s2 = packageVersion("s2") > 
    ##     "1.0.7", k = 200, dwithin = TRUE) 
    ## NULL

If we permit **s2** methods to run, without other arguments set, and
**s2** \> 1.0-7,
[`s2::s2_dwithin_matrix()`](https://r-spatial.github.io/s2/reference/s2_closest_feature.html)
is run:

``` r
if (packageVersion("s2") > "1.0.7") {
  system.time(for (i in 1:(reps/5)) suppressWarnings(pts_ll3_nb <- dnearneigh(pts_ll, d1=0,
      d2=0.75*max_1nn_ll)))/(reps/5)
}
```

    ##    user  system elapsed 
    ##  0.0520  0.0000  0.0525

Alternatively, spherical distances can be used with `dwithin=FALSE` and
[`s2::s2_closest_edges()`](https://r-spatial.github.io/s2/reference/s2_closest_feature.html);
although running in similar time,
[`s2::s2_closest_edges()`](https://r-spatial.github.io/s2/reference/s2_closest_feature.html)
depends on the additional `k=` argument, which, if mis-set, may miss
valid neighbours:

``` r
system.time(for (i in 1:(reps/5)) suppressWarnings(pts_ll5_nb <- dnearneigh(pts_ll, d1=0, d2=0.75*max_1nn_ll, dwithin=FALSE)))/(reps/5)
```

    ##    user  system elapsed 
    ##   0.041   0.000   0.041

``` r
if (packageVersion("s2") > "1.0.7") all.equal(pts_ll3_nb, pts_ll5_nb, check.attributes=FALSE)
```

    ## [1] TRUE

Using
[`s2::s2_closest_edges()`](https://r-spatial.github.io/s2/reference/s2_closest_feature.html)
respects `d1 > 0` without requiring a second pass in R, so is faster
than
[`s2::s2_dwithin_matrix()`](https://r-spatial.github.io/s2/reference/s2_closest_feature.html):

``` r
if (packageVersion("s2") > "1.0.7") {
  system.time(for (i in 1:(reps/5)) suppressWarnings(pts_ll3a_nb <- dnearneigh(pts_ll, d1=5,
      d2=0.75*max_1nn_ll, dwithin=FALSE)))/(reps/5)
}
```

    ##    user  system elapsed 
    ##  0.0380  0.0005  0.0385

Using
[`s2::s2_dwithin_matrix()`](https://r-spatial.github.io/s2/reference/s2_closest_feature.html)
requires a second pass, one for the lower bound, another for the upper
bound, and a set difference operation to find neighbours in the distance
band:

``` r
if (packageVersion("s2") > "1.0.7") {
    system.time(for (i in 1:(reps/5)) suppressWarnings(pts_ll5a_nb <- dnearneigh(pts_ll, d1=5,
        d2=0.75*max_1nn_ll)))/(reps/5)
}
```

    ##    user  system elapsed 
    ##  0.0850  0.0000  0.0855

``` r
if (packageVersion("s2") > "1.0.7") all.equal(pts_ll3a_nb, pts_ll5a_nb, check.attributes=FALSE)
```

    ## [1] TRUE

Setting `use_s2=FALSE` falls back to the legacy version, which uses
symmetry to reduce time:

``` r
system.time(for (i in 1:reps) suppressWarnings(pts_ll6_nb <- dnearneigh(pts_ll, d1=0, d2=0.75*max_1nn_ll, use_s2=FALSE)))/reps
```

    ##    user  system elapsed 
    ##  0.0141  0.0001  0.0143

Minor differences may occur between the legacy ellipsoid and **s2**
spherical approaches:

``` r
all.equal(pts_ll5_nb, pts_ll6_nb, check.attributes=FALSE)
```

    ##  [1] "Component 20: Numeric: lengths (6, 5) differ"     
    ##  [2] "Component 28: Numeric: lengths (7, 6) differ"     
    ##  [3] "Component 112: Numeric: lengths (109, 108) differ"
    ##  [4] "Component 116: Numeric: lengths (109, 108) differ"
    ##  [5] "Component 122: Numeric: lengths (105, 106) differ"
    ##  [6] "Component 123: Numeric: lengths (108, 107) differ"
    ##  [7] "Component 130: Numeric: lengths (108, 109) differ"
    ##  [8] "Component 134: Numeric: lengths (106, 105) differ"
    ##  [9] "Component 158: Numeric: lengths (101, 102) differ"
    ## [10] "Component 165: Numeric: lengths (101, 102) differ"
    ## [11] "Component 168: Numeric: lengths (101, 102) differ"
    ## [12] "Component 179: Numeric: lengths (89, 90) differ"  
    ## [13] "Component 180: Numeric: lengths (96, 97) differ"  
    ## [14] "Component 188: Numeric: lengths (46, 47) differ"  
    ## [15] "Component 189: Numeric: lengths (55, 56) differ"  
    ## [16] "Component 196: Numeric: lengths (47, 46) differ"  
    ## [17] "Component 210: Numeric: lengths (106, 104) differ"
    ## [18] "Component 226: Numeric: lengths (88, 87) differ"  
    ## [19] "Component 229: Numeric: lengths (55, 53) differ"  
    ## [20] "Component 235: Numeric: lengths (40, 39) differ"  
    ## [21] "Component 237: Numeric: lengths (14, 15) differ"  
    ## [22] "Component 245: Numeric: lengths (16, 15) differ"

``` r
system.time(for (i in 1:reps) suppressWarnings(pts_ll6a_nb <- dnearneigh(pts_ll, d1=5, d2=0.75*max_1nn_ll, use_s2=FALSE)))/reps
```

    ##    user  system elapsed 
    ##  0.0138  0.0000  0.0139

``` r
if (packageVersion("s2") > "1.0.7") all.equal(pts_ll5a_nb, pts_ll6a_nb, check.attributes=FALSE)
```

    ##  [1] "Component 20: Numeric: lengths (6, 5) differ"       
    ##  [2] "Component 28: Numeric: lengths (7, 6) differ"       
    ##  [3] "Component 112: Numeric: lengths (62, 61) differ"    
    ##  [4] "Component 113: Numeric: lengths (62, 63) differ"    
    ##  [5] "Component 116: Numeric: lengths (56, 55) differ"    
    ##  [6] "Component 119: Numeric: lengths (68, 69) differ"    
    ##  [7] "Component 122: Numeric: lengths (43, 44) differ"    
    ##  [8] "Component 123: Numeric: lengths (50, 49) differ"    
    ##  [9] "Component 128: Numeric: lengths (65, 64) differ"    
    ## [10] "Component 130: Numeric: lengths (61, 63) differ"    
    ## [11] "Component 132: Numeric: lengths (45, 46) differ"    
    ## [12] "Component 134: Numeric: lengths (46, 45) differ"    
    ## [13] "Component 136: Numeric: lengths (61, 62) differ"    
    ## [14] "Component 147: Numeric: lengths (50, 51) differ"    
    ## [15] "Component 154: Numeric: lengths (56, 57) differ"    
    ## [16] "Component 158: Numeric: lengths (49, 50) differ"    
    ## [17] "Component 165: Mean relative difference: 0.02823018"
    ## [18] "Component 168: Numeric: lengths (54, 56) differ"    
    ## [19] "Component 179: Numeric: lengths (77, 78) differ"    
    ## [20] "Component 180: Numeric: lengths (85, 86) differ"    
    ## [21] "Component 188: Numeric: lengths (39, 40) differ"    
    ## [22] "Component 189: Numeric: lengths (48, 49) differ"    
    ## [23] "Component 196: Numeric: lengths (45, 44) differ"    
    ## [24] "Component 210: Numeric: lengths (68, 66) differ"    
    ## [25] "Component 226: Numeric: lengths (82, 81) differ"    
    ## [26] "Component 229: Numeric: lengths (48, 46) differ"    
    ## [27] "Component 235: Numeric: lengths (38, 37) differ"    
    ## [28] "Component 237: Numeric: lengths (14, 15) differ"    
    ## [29] "Component 245: Numeric: lengths (15, 14) differ"

#### Contiguity neighbours for spherical polygon support

It also turns out that when **sf** functions are used to find contiguity
neighbours, **s2** spatial indexing functionality is accessed in finding
candidate neighbours in intersecting geometries.

``` r
NY8_sf_ll <- st_transform(NY8_sf, "OGC:CRS84")
st_is_longlat(NY8_sf_ll)
```

    ## [1] TRUE

The timings are a little slower when
[`st_intersects()`](https://r-spatial.github.io/sf/reference/geos_binary_pred.html)
hands off geometry predicates to `s2_intersects_matrix()`, but the
results are the same, and because spatial indexing is used, this scales
well for larger data sets:

``` r
sf_use_s2(TRUE)
system.time(for (i in 1:reps) NY8_sf_1_nb_ll <- poly2nb(NY8_sf_ll, queen=TRUE, snap=eps))/reps
```

    ##    user  system elapsed 
    ##  0.1659  0.0024  0.1689

``` r
all.equal(NY8_sf_1_nb, NY8_sf_1_nb_ll, check.attributes=FALSE)
```

    ## [1] TRUE
