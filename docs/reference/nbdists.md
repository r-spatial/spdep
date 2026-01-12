# Spatial link distance measures

Given a list of spatial neighbour links (a neighbours list of object
type `nb`), the function returns the Euclidean distances along the links
in a list of the same form as the neighbours list. If longlat = TRUE,
Great Circle distances are used.

## Usage

``` r
nbdists(nb, coords, longlat = NULL)
```

## Arguments

- nb:

  an object of class `nb`

- coords:

  matrix of point coordinates, an object inheriting from SpatialPoints
  or an `"sf"` or `"sfc"` object; if the `"sf"` or `"sfc"` object
  geometries are in geographical coordinates
  (`sf::st_is_longlat(x) == TRUE` and `sf::sf_use_s2() == TRUE`), s2
  will be used to find distances
  <https://github.com/r-spatial/s2/issues/125>

- longlat:

  TRUE if point coordinates are longitude-latitude decimal degrees, in
  which case distances are measured in kilometers; if coords is a
  SpatialPoints object, the value is taken from the object itself

## Value

A list with class `nbdist`

## Author

Roger Bivand <Roger.Bivand@nhh.no>

## See also

[`summary.nb`](https://r-spatial.github.io/spdep/reference/summary.nb.md),
[`nb2listw`](https://r-spatial.github.io/spdep/reference/nb2listw.md)

## Examples

``` r
columbus <- st_read(system.file("shapes/columbus.gpkg", package="spData")[1], quiet=TRUE)
col.gal.nb <- read.gal(system.file("weights/columbus.gal", package="spData")[1])
coords <- st_coordinates(st_centroid(columbus))
#> Warning: st_centroid assumes attributes are constant over geometries
dlist <- nbdists(col.gal.nb, coords)
dlist <- lapply(dlist, function(x) 1/x)
stem(unlist(dlist))
#> 
#>   The decimal point is at the |
#> 
#>   1 | 1111112222223333334444444444
#>   1 | 555555556666667777777777888888888899999999999999999999
#>   2 | 00000000111111111111111122222222222222333333333333334444
#>   2 | 55555555556666666666666677777777888888999999
#>   3 | 0000111111222233334444
#>   3 | 77779999
#>   4 | 223344
#>   4 | 77
#>   5 | 001144
#>   5 | 
#>   6 | 
#>   6 | 99
#>   7 | 
#>   7 | 88
#> 
```
