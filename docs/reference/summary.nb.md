# Print and summary function for neighbours and weights lists

The function prints summary measures for links in a neighbours list. If
a matrix of coordinates is given as well, summary descriptive measures
for the link lengths are also printed. Print and summary functions are
also available for `"listw"` weights list objects, also reporting
constants (S0, S1, S2) used in inference for global spatial
autocorrelation statistics such as Moran's I, Geary's C, join-count
tests and Getis-Ord G.

## Usage

``` r
# S3 method for class 'nb'
summary(object, coords=NULL, longlat = NULL, scale = 1, ...)
# S3 method for class 'nb'
print(x, ...)
# S3 method for class 'listw'
summary(object, coords, longlat, zero.policy = attr(object, "zero.policy"), 
 scale = 1, adjust.n=TRUE, ...)
# S3 method for class 'listw'
print(x, zero.policy = attr(x, "zero.policy"), ...)
```

## Arguments

- object:

  an object of class `nb`

- coords:

  matrix of region point coordinates or a SpatialPoints object or an
  `sfc` points object

- longlat:

  TRUE if point coordinates are longitude-latitude decimal degrees, in
  which case distances are measured in kilometers; if coords is a
  SpatialPoints object, the value is taken from the object itself

- ...:

  additional arguments affecting the output produced

- x:

  an object of class `nb`

- zero.policy:

  default `attr(listw, "zero.policy")` as set when `listw` was created,
  if attribute not set, use global option value; if FALSE stop with
  error for any empty neighbour sets

- scale:

  passed through to [`stem()`](https://rdrr.io/r/graphics/stem.html) for
  control of plot length

- adjust.n:

  default TRUE, if FALSE the number of observations is not adjusted for
  no-neighbour observations, if TRUE, the number of observations in
  `spweights.constants` is adjusted

## Author

Roger Bivand <Roger.Bivand@nhh.no>

## See also

[`plot.nb`](https://r-spatial.github.io/spdep/reference/plot.nb.md)

## Examples

``` r
columbus <- st_read(system.file("shapes/columbus.gpkg", package="spData")[1], quiet=TRUE)
col.gal.nb <- read.gal(system.file("weights/columbus.gal", package="spData")[1])
coords <- st_centroid(st_geometry(columbus), of_largest_polygon=TRUE)
col.gal.nb
#> Neighbour list object:
#> Number of regions: 49 
#> Number of nonzero links: 230 
#> Percentage nonzero weights: 9.579342 
#> Average number of links: 4.693878 
summary(col.gal.nb, coords)
#> Neighbour list object:
#> Number of regions: 49 
#> Number of nonzero links: 230 
#> Percentage nonzero weights: 9.579342 
#> Average number of links: 4.693878 
#> Link number distribution:
#> 
#>  2  3  4  5  6  7  8  9 10 
#>  7  7 13  4  9  6  1  1  1 
#> 7 least connected regions:
#> 1 6 31 39 42 46 47 with 2 links
#> 1 most connected region:
#> 20 with 10 links
col.listw <- nb2listw(col.gal.nb, style="W")
col.listw
#> Characteristics of weights list object:
#> Neighbour list object:
#> Number of regions: 49 
#> Number of nonzero links: 230 
#> Percentage nonzero weights: 9.579342 
#> Average number of links: 4.693878 
#> 
#> Weights style: W 
#> Weights constants summary:
#>    n   nn S0       S1       S2
#> W 49 2401 49 23.48489 204.6687
summary(col.listw)
#> Characteristics of weights list object:
#> Neighbour list object:
#> Number of regions: 49 
#> Number of nonzero links: 230 
#> Percentage nonzero weights: 9.579342 
#> Average number of links: 4.693878 
#> Link number distribution:
#> 
#>  2  3  4  5  6  7  8  9 10 
#>  7  7 13  4  9  6  1  1  1 
#> 7 least connected regions:
#> 1 6 31 39 42 46 47 with 2 links
#> 1 most connected region:
#> 20 with 10 links
#> 
#> Weights style: W 
#> Weights constants summary:
#>    n   nn S0       S1       S2
#> W 49 2401 49 23.48489 204.6687
col_geoms <- st_geometry(columbus)
col_geoms[21] <- st_buffer(col_geoms[21], dist=-0.05)
st_geometry(columbus) <- col_geoms
nb <- poly2nb(columbus)
#> Warning: some observations have no neighbours;
#> if this seems unexpected, try increasing the snap argument.
#> Warning: neighbour object has 3 sub-graphs;
#> if this sub-graph count seems unexpected, try increasing the snap argument.
summary(nb)
#> Neighbour list object:
#> Number of regions: 49 
#> Number of nonzero links: 230 
#> Percentage nonzero weights: 9.579342 
#> Average number of links: 4.693878 
#> 1 region with no links:
#> 21
#> 3 disjoint connected subgraphs
#> Link number distribution:
#> 
#>  0  2  3  4  5  6  7  8  9 10 
#>  1  5  9 12  4 10  2  4  1  1 
#> 5 least connected regions:
#> 1 6 42 46 47 with 2 links
#> 1 most connected region:
#> 20 with 10 links
try(nb2listw(nb, style="W"))
#> Error in nb2listw(nb, style = "W") : 
#>   Empty neighbour sets found (zero.policy: FALSE)
summary(nb2listw(nb, style="W", zero.policy=TRUE))
#> Characteristics of weights list object:
#> Neighbour list object:
#> Number of regions: 49 
#> Number of nonzero links: 230 
#> Percentage nonzero weights: 9.579342 
#> Average number of links: 4.693878 
#> 1 region with no links:
#> 21
#> 3 disjoint connected subgraphs
#> Link number distribution:
#> 
#>  0  2  3  4  5  6  7  8  9 10 
#>  1  5  9 12  4 10  2  4  1  1 
#> 5 least connected regions:
#> 1 6 42 46 47 with 2 links
#> 1 most connected region:
#> 20 with 10 links
#> 
#> Weights style: W 
#> Weights constants summary:
#>    n   nn S0       S1       S2
#> W 48 2304 48 22.46811 199.4398
summary(nb2listw(nb, style="W", zero.policy=TRUE), adjust.n=FALSE)
#> Characteristics of weights list object:
#> Neighbour list object:
#> Number of regions: 49 
#> Number of nonzero links: 230 
#> Percentage nonzero weights: 9.579342 
#> Average number of links: 4.693878 
#> 1 region with no links:
#> 21
#> 3 disjoint connected subgraphs
#> Link number distribution:
#> 
#>  0  2  3  4  5  6  7  8  9 10 
#>  1  5  9 12  4 10  2  4  1  1 
#> 5 least connected regions:
#> 1 6 42 46 47 with 2 links
#> 1 most connected region:
#> 20 with 10 links
#> 
#> Weights style: W 
#> Weights constants summary:
#>    n   nn S0       S1       S2
#> W 49 2401 48 22.46811 199.4398
```
