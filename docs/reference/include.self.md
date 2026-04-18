# Include self in neighbours list

The function includes the region itself in its own list of neighbours,
and sets attribute "self.included" to TRUE; `remove.self` reverts the
effects of `include.self`.

## Usage

``` r
include.self(nb)
remove.self(nb)
```

## Arguments

- nb:

  input neighbours list of class `nb`

## Value

The function returns an object of class `nb` with a list of integer
vectors containing neighbour region number ids; attribute
"self.included" is set to TRUE.

## Author

Roger Bivand <Roger.Bivand@nhh.no>

## See also

[`summary.nb`](https://r-spatial.github.io/spdep/reference/summary.nb.md)

## Examples

``` r
columbus <- st_read(system.file("shapes/columbus.gpkg", package="spData")[1], quiet=TRUE)
col.gal.nb <- read.gal(system.file("weights/columbus.gal", package="spData")[1])
coords <- st_coordinates(st_centroid(columbus))
#> Warning: st_centroid assumes attributes are constant over geometries
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
summary(include.self(col.gal.nb), coords)
#> Neighbour list object:
#> Number of regions: 49 
#> Number of nonzero links: 279 
#> Percentage nonzero weights: 11.62016 
#> Average number of links: 5.693878 
#> Link number distribution:
#> 
#>  3  4  5  6  7  8  9 10 11 
#>  7  7 13  4  9  6  1  1  1 
#> 7 least connected regions:
#> 1 6 31 39 42 46 47 with 3 links
#> 1 most connected region:
#> 20 with 11 links
summary(remove.self(include.self(col.gal.nb)), coords)
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
```
