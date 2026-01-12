# Convert a square spatial weights matrix to a weights list object

The function converts a square spatial weights matrix, optionally a
sparse matrix to a weights list object, optionally adding region IDs
from the row names of the matrix, as a sequence of numbers 1:nrow(x), or
as given as an argument. The style can be imposed by rebuilting the
weights list object internally.

## Usage

``` r
mat2listw(x, row.names = NULL, style=NULL, zero.policy = NULL)
```

## Arguments

- x:

  A square non-negative matrix with no NAs representing spatial weights;
  may be a matrix of class “sparseMatrix”

- row.names:

  row names to use for region IDs

- style:

  default NULL, missing, set to "M" and warning given; if not "M",
  passed to
  [`nb2listw`](https://r-spatial.github.io/spdep/reference/nb2listw.md)
  to re-build the object

- zero.policy:

  default NULL, use global option value; if FALSE stop with error for
  any empty neighbour sets, if TRUE permit the weights list to be formed
  with zero-length weights vectors

## Value

A `listw` object with the following members:

- style:

  "M", meaning matrix style, underlying style unknown, or assigned style
  argument in rebuilt object

- neighbours:

  the derived neighbours list

- weights:

  the weights for the neighbours derived from the matrix

## Author

Roger Bivand <Roger.Bivand@nhh.no>

## See also

[`nb2listw`](https://r-spatial.github.io/spdep/reference/nb2listw.md),
[`nb2mat`](https://r-spatial.github.io/spdep/reference/nb2mat.md)

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
col005.w.mat <- nb2mat(col005, style="W", zero.policy=TRUE)
try(col005.w.b <- mat2listw(col005.w.mat, style="W"))
#> Warning: no-neighbour observations found, set zero.policy to TRUE;
#> this warning will soon become an error
#> Warning: neighbour object has 8 sub-graphs
#> Error in nb2listw(res$neighbours, glist = res$weights, style = style,  : 
#>   Empty neighbour sets found (zero.policy: FALSE)
col005.w.b <- mat2listw(col005.w.mat, style="W", zero.policy=TRUE)
#> Warning: neighbour object has 8 sub-graphs
summary(col005.w.b$neighbours)
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
diffnb(col005, col005.w.b$neighbours)
#> Warning: neighbour object has 49 sub-graphs
#> Neighbour list object:
#> Number of regions: 49 
#> Number of nonzero links: 0 
#> Percentage nonzero weights: 0 
#> Average number of links: 0 
#> 49 regions with no links:
#> 1005, 1001, 1006, 1002, 1007, 1008, 1004, 1003, 1018, 1010, 1038, 1037,
#> 1039, 1040, 1009, 1036, 1011, 1042, 1041, 1017, 1043, 1019, 1012, 1035,
#> 1032, 1020, 1021, 1031, 1033, 1034, 1045, 1013, 1022, 1044, 1023, 1046,
#> 1030, 1024, 1047, 1016, 1014, 1049, 1029, 1025, 1028, 1048, 1015, 1027,
#> 1026
#> 49 disjoint connected subgraphs
col005.w.mat.3T <- kronecker(diag(3), col005.w.mat)
col005.w.b.3T <- mat2listw(col005.w.mat.3T, style="W", zero.policy=TRUE)
#> Warning: neighbour object has 24 sub-graphs
summary(col005.w.b.3T$neighbours)
#> Neighbour list object:
#> Number of regions: 147 
#> Number of nonzero links: 510 
#> Percentage nonzero weights: 2.360128 
#> Average number of links: 3.469388 
#> 12 regions with no links:
#> 1, 3, 6, 21, 50, 52, 55, 70, 99, 101, 104, 119
#> 24 disjoint connected subgraphs
#> Link number distribution:
#> 
#>  0  1  2  3  4  5  6  7  8  9 
#> 12 33 15 24  9 27  6  6  9  6 
#> 33 least connected regions:
#> 2 5 9 10 31 34 36 39 42 46 47 51 54 58 59 80 83 85 88 91 95 96 100 103 107 108 129 132 134 137 140 144 145 with 1 link
#> 6 most connected regions:
#> 11 16 60 65 109 114 with 9 links
run <- FALSE
if (require("spatialreg", quiet=TRUE)) run <- TRUE
if (run) {
W <- as(nb2listw(col005, style="W", zero.policy=TRUE), "CsparseMatrix")
try(col005.spM <- mat2listw(W))
col005.spM <- mat2listw(W, style="W", zero.policy=TRUE)
summary(col005.spM$neighbours)
}
#> Warning: style is M (missing); style should be set to a valid value
#> Warning: style is M (missing); style should be set to a valid value
#> Warning: 1005, 1006, 1008, 1043 are not origins
#> Warning: neighbour object has 8 sub-graphs
#> Warning: no-neighbour observations found, set zero.policy to TRUE;
#> this warning will soon become an error
#> Warning: neighbour object has 8 sub-graphs
#> Warning: neighbour object has 8 sub-graphs
#> Warning: neighbour object has 8 sub-graphs
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
if (run) {
diffnb(col005, col005.spM$neighbours)
}
#> Warning: neighbour object has 49 sub-graphs
#> Neighbour list object:
#> Number of regions: 49 
#> Number of nonzero links: 0 
#> Percentage nonzero weights: 0 
#> Average number of links: 0 
#> 49 regions with no links:
#> 1005, 1001, 1006, 1002, 1007, 1008, 1004, 1003, 1018, 1010, 1038, 1037,
#> 1039, 1040, 1009, 1036, 1011, 1042, 1041, 1017, 1043, 1019, 1012, 1035,
#> 1032, 1020, 1021, 1031, 1033, 1034, 1045, 1013, 1022, 1044, 1023, 1046,
#> 1030, 1024, 1047, 1016, 1014, 1049, 1029, 1025, 1028, 1048, 1015, 1027,
#> 1026
#> 49 disjoint connected subgraphs
if (run && require("Matrix", quiet=TRUE)) {
IW <- kronecker(Diagonal(3), W)
col005.spM.3T <- mat2listw(as(IW, "CsparseMatrix"), style="W", zero.policy=TRUE)
summary(col005.spM.3T$neighbours)
}
#> Warning: neighbour object has 24 sub-graphs
#> Warning: neighbour object has 24 sub-graphs
#> Neighbour list object:
#> Number of regions: 147 
#> Number of nonzero links: 510 
#> Percentage nonzero weights: 2.360128 
#> Average number of links: 3.469388 
#> 12 regions with no links:
#> 1, 3, 6, 21, 50, 52, 55, 70, 99, 101, 104, 119
#> 24 disjoint connected subgraphs
#> Link number distribution:
#> 
#>  0  1  2  3  4  5  6  7  8  9 
#> 12 33 15 24  9 27  6  6  9  6 
#> 33 least connected regions:
#> 2 5 9 10 31 34 36 39 42 46 47 51 54 58 59 80 83 85 88 91 95 96 100 103 107 108 129 132 134 137 140 144 145 with 1 link
#> 6 most connected regions:
#> 11 16 60 65 109 114 with 9 links
```
