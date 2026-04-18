# Subset a neighbours list

The function subsets a neighbors list, retaining objects for which the
subset argument vector is TRUE.

## Usage

``` r
# S3 method for class 'nb'
subset(x, subset, ...)
```

## Arguments

- x:

  an object of class `nb`

- subset:

  logical expression

- ...:

  generic function pass-through

## Value

The function returns an object of class `nb` with a list of integer
vectors containing neighbour region number ids (compacted to run from
1:number of regions in subset).

## Author

Roger Bivand <Roger.Bivand@nhh.no>

## See also

[`nb2listw`](https://r-spatial.github.io/spdep/reference/nb2listw.md)

## Examples

``` r
columbus <- st_read(system.file("shapes/columbus.gpkg", package="spData")[1], quiet=TRUE)
col.gal.nb <- read.gal(system.file("weights/columbus.gal", package="spData")[1])
coords <- st_coordinates(st_centroid(columbus))
#> Warning: st_centroid assumes attributes are constant over geometries
plot(col.gal.nb, coords)
to.be.dropped <- c(31, 34, 36, 39, 42, 46)
text(coords[to.be.dropped,1], coords[to.be.dropped,2], labels=to.be.dropped,
  pos=2, offset=0.3)
sub.col.gal.nb <- subset(col.gal.nb,
  !(1:length(col.gal.nb) %in% to.be.dropped))
plot(sub.col.gal.nb, coords[-to.be.dropped,], col="red", add=TRUE)

which(!(attr(col.gal.nb, "region.id") %in%
  attr(sub.col.gal.nb, "region.id")))
#> [1] 31 34 36 39 42 46
```
