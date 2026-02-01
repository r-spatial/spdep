# Data set with 4 life condition indices of Belo Horizonte region

The data are collected inthe Atlas of condition indices published by the
Joao Pinheiro Foundation and UNDP.

## Format

A shape polygon object with seven variables:

- id:

  The identificator

- Name:

  Name of city

- Population:

  The population of city

- HLCI:

  Health Life Condition Index

- ELCI:

  Education Life Condition Index

- CLCI:

  Children Life Condition Index

- ELCI:

  Economic Life Condition Index

## Examples

``` r
GDAL37 <- numeric_version(unname(sf::sf_extSoftVersion()["GDAL"]), strict=FALSE)
(GDAL37 <- ifelse(is.na(GDAL37), FALSE, GDAL37 >= "3.7.0"))
#> [1] TRUE
file <- "etc/shapes/bhicv.gpkg.zip"
zipfile <- system.file(file, package="spdep")
if (GDAL37) {
    bh <- st_read(zipfile)
} else {
    td <- tempdir()
    bn <- sub(".zip", "", basename(file), fixed=TRUE)
    target <- unzip(zipfile, files=bn, exdir=td)
    bh <- st_read(target)
}
#> Reading layer `bhicv' from data source 
#>   `/tmp/RtmpEjBomy/temp_libpath486f623f1b6ee/spdep/etc/shapes/bhicv.gpkg.zip' 
#>   using driver `GPKG'
#> Simple feature collection with 98 features and 8 fields
#> Geometry type: POLYGON
#> Dimension:     XY
#> Bounding box:  xmin: -45.02175 ymin: -20.93007 xmax: -42.50321 ymax: -18.08342
#> Geodetic CRS:  Corrego Alegre 1970-72
```
