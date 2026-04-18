# Differences between neighbours lists

The function finds symmetric differences between lists of neighbours
(ordering of objects does not matter), returning a `nb` neighbour list
of those found

## Usage

``` r
diffnb(x, y, verbose=NULL, legacy=TRUE)
```

## Arguments

- x:

  an object of class `nb`

- y:

  an object of class `nb`

- verbose:

  default NULL, use global option value; report regions ids taken from
  object attribute "region.id" with differences, ignored when `legacy=`
  is FALSE

- legacy:

  default TRUE; use legacy code; if FALSE use differences between sparse
  matrix representations

## Value

A neighbours list with class `nb`

## Author

Roger Bivand <Roger.Bivand@nhh.no>

## Examples

``` r
columbus <- st_read(system.file("shapes/columbus.gpkg", package="spData")[1], quiet=TRUE)
coords <- st_centroid(st_geometry(columbus), of_largest_polygon=TRUE)
rn <- row.names(columbus)
knn1 <- knearneigh(coords, 1)
knn2 <- knearneigh(coords, 2)
nb1 <- knn2nb(knn1, row.names=rn)
#> Warning: neighbour object has 13 sub-graphs
nb2 <- knn2nb(knn2, row.names=rn)
diffs <- diffnb(nb2, nb1, legacy=TRUE)
#> Warning: neighbour object has 10 sub-graphs
diffsl <- diffnb(nb2, nb1, legacy=FALSE)
#> Warning: neighbour object has 10 sub-graphs
#> Warning: neighbour object has 10 sub-graphs
#> Warning: neighbour object has 10 sub-graphs
all.equal(diffs, diffsl)
#> [1] "Attributes: < Component “call”: target, current do not match when deparsed >"
# call attribute fifference expected
diffsd <- union.nb(setdiff.nb(nb1, nb2), setdiff.nb(nb2, nb1))
#> Warning: neighbour object has 49 sub-graphs
#> Warning: neighbour object has 10 sub-graphs
#> Warning: neighbour object has 10 sub-graphs
all.equal(diffs, diffsd)
#> [1] "Attributes: < Component “call”: target, current do not match when deparsed >"
# call attribute fifference expected
opar <- par(no.readonly=TRUE)
plot(st_geometry(columbus), border="grey", reset=FALSE,
 main="Plot of first (black) and second (red)\nnearest neighbours")
plot(nb1, coords, add=TRUE)
plot(diffs, coords, add=TRUE, col="red", lty=2)

par(opar)
```
