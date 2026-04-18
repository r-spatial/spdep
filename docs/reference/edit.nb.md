# Interactive editing of neighbours lists

The function provides simple interactive editing of neighbours lists to
allow unneeded links to be deleted, and missing links to be inserted. It
uses `identify` to pick the endpoints of the link to be deleted or
added, and asks for confirmation before committing. If the result is not
assigned to a new object, the editing will be lost - as in `edit`.

This method relies on direct contact with the graphics device. Do not
use in RStudio.

## Usage

``` r
# S3 method for class 'nb'
edit(name, coords, polys=NULL, ..., use_region.id=FALSE)
```

## Arguments

- name:

  an object of class `nb`

- coords:

  matrix of region point coordinates; if missing and polys= inherits
  from `SpatialPolygons`, the label points of that object are used

- polys:

  if polygon boundaries supplied, will be used as background; must
  inherit from `SpatialPolygons`

- ...:

  further arguments passed to or from other methods

- use_region.id:

  default `FALSE`, in `identify` use 1-based observation numbers,
  otherwise use the `nb` `region.id` attribute values

## Value

The function returns an object of class `nb` with the edited list of
integer vectors containing neighbour region number ids, with added
attributes tallying the added and deleted links.

## Author

Roger Bivand <Roger.Bivand@nhh.no>

## See also

[`summary.nb`](https://r-spatial.github.io/spdep/reference/summary.nb.md),
[`plot.nb`](https://r-spatial.github.io/spdep/reference/plot.nb.md)

## Examples

``` r
# \dontrun{
columbus <- st_read(system.file("shapes/columbus.gpkg", package="spData")[1], quiet=TRUE)
class(columbus)
#> [1] "sf"         "data.frame"
if (FALSE) nnb1 <- edit.nb(col.gal.nb, polys=as(columbus, "Spatial"))
# }
```
