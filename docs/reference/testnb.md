# Test a neighbours list for symmetry

Checks a neighbours list for symmetry/transitivity (if i is a neighbour
of j, then j is a neighbour of i). This holds for distance and
contiguity based neighbours, but not for k-nearest neighbours. The
helper function `sym.attr.nb()` calls `is.symmetric.nb()` to set the
`sym` attribute if needed, and `make.sym.nb` makes a non-symmetric list
symmetric by adding neighbors. `is.symmetric.glist` checks a list of
general weights corresponding to neighbours for symmetry for symmetric
neighbours.

## Usage

``` r
is.symmetric.nb(nb, verbose = NULL, force = FALSE)
sym.attr.nb(nb)
make.sym.nb(nb)
old.make.sym.nb(nb)
is.symmetric.glist(nb, glist)
```

## Arguments

- nb:

  an object of class `nb` with a list of integer vectors containing
  neighbour region number ids.

- verbose:

  default NULL, use global option value; if TRUE prints non-matching
  pairs

- force:

  do not respect a neighbours list `sym` attribute and test anyway

- glist:

  list of general weights corresponding to neighbours

## Value

TRUE if symmetric, FALSE if not; is.symmetric.glist returns a value with
an attribute, "d", indicating for failed symmetry the largest failing
value.

## Note

A new version of `make.sym.nb` by Bjarke Christensen is now included.
The older version has been renamed `old.make.sym.nb`, and their
comparison constitutes a nice demonstration of vectorising speedup using
`sapply` and `lapply` rather than loops. When any no-neighbour
observations are present, `old.make.sym.nb` is used.

## Author

Roger Bivand <Roger.Bivand@nhh.no>

## See also

[`read.gal`](https://r-spatial.github.io/spdep/reference/read.gal.md)

## Examples

``` r
columbus <- st_read(system.file("shapes/columbus.gpkg", package="spData")[1], quiet=TRUE)
col.gal.nb <- read.gal(system.file("weights/columbus.gal", package="spData")[1])
coords <- st_coordinates(st_centroid(columbus))
#> Warning: st_centroid assumes attributes are constant over geometries
ind <- row.names(as(columbus, "Spatial"))
print(is.symmetric.nb(col.gal.nb, verbose=TRUE, force=TRUE))
#> [1] TRUE
k4 <- knn2nb(knearneigh(coords, k=4), row.names=ind)
k4 <- sym.attr.nb(k4)
print(is.symmetric.nb(k4))
#> [1] FALSE
k4.sym <- make.sym.nb(k4)
print(is.symmetric.nb(k4.sym))
#> [1] TRUE
```
