# Columbus OH spatial analysis data set

The data set is now part of the spData package

## Usage

``` r
data(columbus)
```

## Examples

``` r
columbus <- st_read(system.file("shapes/columbus.gpkg", package="spData")[1], quiet=TRUE)
col.gal.nb <- read.gal(system.file("weights/columbus.gal", package="spData")[1])
```
