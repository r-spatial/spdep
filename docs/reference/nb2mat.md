# Spatial weights matrices for neighbours lists

The function generates a weights matrix for a neighbours list with
spatial weights for the chosen coding scheme.

## Usage

``` r
nb2mat(neighbours, glist=NULL, style="W", zero.policy=NULL)
listw2mat(listw)
```

## Arguments

- neighbours:

  an object of class `nb`

- glist:

  list of general weights corresponding to neighbours

- style:

  `style` can take values W, B, C, and S

- zero.policy:

  default NULL, use global option value; if FALSE stop with error for
  any empty neighbour sets, if TRUE permit the weights list to be formed
  with zero-length weights vectors

- listw:

  a `listw` object from for example `nb2listw`

## Details

Starting from a binary neighbours list, in which regions are either
listed as neighbours or are absent (thus not in the set of neighbours
for some definition), the function creates an n by n weights matrix with
values given by the coding scheme style chosen. B is the basic binary
coding, W is row standardised, C is globally standardised, while S is
the variance-stabilizing coding scheme proposed by Tiefelsdorf et al.
1999, p. 167-168.

The function leaves matrix rows as zero for any regions with zero
neighbours fore zero.policy TRUE. These will in turn generate lag values
of zero, equivalent to the sum of products of the zero row
`t(rep(0, length=length(neighbours))) %*% x`, for arbitraty numerical
vector `x` of length `length(neighbours)`. The spatially lagged value of
x for the zero-neighbour region will then be zero, which may (or may
not) be a sensible choice.

## Value

An n by n matrix, where n=length(neighbours)

## References

Tiefelsdorf, M., Griffith, D. A., Boots, B. 1999 A variance-stabilizing
coding scheme for spatial link matrices, Environment and Planning A, 31,
pp. 165-180.

## Author

Roger Bivand <Roger.Bivand@nhh.no>

## See also

[`nb2listw`](https://r-spatial.github.io/spdep/reference/nb2listw.md)

## Examples

``` r
columbus <- st_read(system.file("shapes/columbus.gpkg", package="spData")[1], quiet=TRUE)
col005 <- dnearneigh(st_coordinates(st_centroid(st_geometry(columbus),
 of_largest_polygon=TRUE)), 0, 0.5, as.character(columbus$NEIGNO))
#> Warning: neighbour object has 8 sub-graphs
summary(col005)
#> Neighbour list object:
#> Number of regions: 49 
#> Number of nonzero links: 170 
#> Percentage nonzero weights: 7.080383 
#> Average number of links: 3.469388 
#> 4 regions with no links:
#> 1005, 1006, 1008, 1043
#> 8 disjoint connected subgraphs
#> Link number distribution:
#> 
#>  0  1  2  3  4  5  6  7  8  9 
#>  4 11  5  8  3  9  2  2  3  2 
#> 11 least connected regions:
#> 1001 1007 1018 1010 1045 1044 1046 1047 1049 1048 1015 with 1 link
#> 2 most connected regions:
#> 1038 1036 with 9 links
col005.w.mat <- nb2mat(col005, style="B", zero.policy=TRUE)
table(round(rowSums(col005.w.mat)))
#> 
#>  0  1  2  3  4  5  6  7  8  9 
#>  4 11  5  8  3  9  2  2  3  2 
```
