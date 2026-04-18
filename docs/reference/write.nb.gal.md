# Write a neighbours list as a GAL lattice file

Write a neighbours list as a GAL lattice file, may also use newer GeoDa
header format

## Usage

``` r
write.nb.gal(nb, file, oldstyle=TRUE, shpfile=NULL, ind=NULL)
```

## Arguments

- nb:

  an object of class `nb` with a list of integer vectors containing
  neighbour region number ids.

- file:

  name of file with GAL lattice data

- oldstyle:

  if TRUE, first line of file contains only number of spatial units, if
  FALSE, uses newer GeoDa style

- shpfile:

  Shapefile name taken from GAL file for this dataset

- ind:

  region id indicator variable name

## Author

Roger Bivand <Roger.Bivand@nhh.no>

## See also

[`read.gal`](https://r-spatial.github.io/spdep/reference/read.gal.md)

## Examples

``` r
col.gal.nb <- read.gal(system.file("weights/columbus.gal", package="spData")[1])
GALfile <- tempfile("GAL")
write.nb.gal(col.gal.nb, GALfile)
col.queen <- read.gal(GALfile)
summary(diffnb(col.queen, col.gal.nb))
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
#> Link number distribution:
#> 
#>  0 
#> 49 
```
