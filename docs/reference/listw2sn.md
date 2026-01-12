# Spatial neighbour sparse representation

The function makes a `"spatial neighbour"` object representation
(similar to the S-PLUS spatial statististics module representation of a
`"listw"` spatial weights object. `sn2listw()` is the inverse function
to `listw2sn()`, creating a `"listw"` object from a
`"spatial neighbour"` object.

## Usage

``` r
listw2sn(listw)
sn2listw(sn, style = NULL, zero.policy = NULL, from_mat2listw=FALSE)
```

## Arguments

- listw:

  a `listw` object from for example `nb2listw`

- sn:

  a `spatial.neighbour` object

- style:

  default NULL, missing, set to "M" and warning given; if not "M",
  passed to
  [`nb2listw`](https://r-spatial.github.io/spdep/reference/nb2listw.md)
  to re-build the object

- zero.policy:

  default NULL, use global option value; if FALSE stop with error for
  any empty neighbour sets, if TRUE permit the weights list to be formed
  with zero-length weights vectors

- from_mat2listw:

  default FALSE, set TRUE if called from `mat2listw`

## Value

`listw2sn()`returns a data frame with three columns, and with class
`spatial.neighbour`:

- from:

  region number id for the start of the link (S-PLUS row.id)

- to:

  region number id for the end of the link (S-PLUS col.id)

- weights:

  weight for this link

## Author

Roger Bivand <Roger.Bivand@nhh.no>

## See also

[`nb2listw`](https://r-spatial.github.io/spdep/reference/nb2listw.md)

## Examples

``` r
columbus <- st_read(system.file("shapes/columbus.gpkg", package="spData")[1], quiet=TRUE)
col.gal.nb <- read.gal(system.file("weights/columbus.gal", package="spData")[1])
col.listw <- nb2listw(col.gal.nb)
col.listw$neighbours[[1]]
#> [1] 2 3
col.listw$weights[[1]]
#> [1] 0.5 0.5
col.sn <- listw2sn(col.listw)
str(col.sn)
#> Classes ‘spatial.neighbour’ and 'data.frame':    230 obs. of  3 variables:
#>  $ from   : int  1 1 2 2 2 3 3 3 3 4 ...
#>  $ to     : int  2 3 1 3 4 1 2 4 5 2 ...
#>  $ weights: num  0.5 0.5 0.333 0.333 0.333 ...
#>  - attr(*, "n")= int 49
#>  - attr(*, "region.id")= chr [1:49] "1" "2" "3" "4" ...
#>  - attr(*, "neighbours.attrs")= chr [1:7] "class" "region.id" "GeoDa" "gal" ...
#>  - attr(*, "weights.attrs")= chr [1:3] "mode" "W" "comp"
#>  - attr(*, "GeoDa")=List of 2
#>   ..$ shpfile: chr NA
#>   ..$ ind    : chr NA
#>  - attr(*, "listw.call")= language nb2listw(neighbours = col.gal.nb)
```
