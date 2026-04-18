# Use vector files for import and export of weights

Use vector files for import and export of weights, storing spatial
entity coordinates in the arcs, and the entity indices in the data
frame.

## Usage

``` r
nb2lines(nb, wts, coords, proj4string=NULL, as_sf=FALSE)
listw2lines(listw, coords, proj4string=NULL, as_sf=FALSE)
df2sn(df, i="i", i_ID="i_ID", j="j", wt="wt")
```

## Arguments

- nb:

  a neighbour object of class `nb`

- wts:

  list of general weights corresponding to neighbours

- coords:

  matrix of region point coordinates, a `Spatial` object (points or
  polygons), or an `sfc` object (points or polygons)

- proj4string:

  default NULL; if `coords` is a Spatial or sf object, this value will
  be used, otherwise the value will be converted appropriately

- as_sf:

  output object in `Spatial` or `sf` format, default FALSE, set to TRUE
  if coords is an `sfc` object and FALSE if a `Spatial` object

- listw:

  a `listw` object of spatial weights

- df:

  a data frame read from a shapefile, derived from the output of
  `nb2lines`

- i:

  character name of column in df with from entity index

- i_ID:

  character name of column in df with from entity region ID

- j:

  character name of column in df with to entity index

- wt:

  character name of column in df with weights

## Details

The neighbour and weights objects may be retrieved by converting the
specified columns of the data slot of the SpatialLinesDataFrame object
into a spatial.neighbour object, which is then converted into a weights
list object.

## Value

`nb2lines` and `listw2lines` return a SpatialLinesDataFrame object or an
sf object; the data frame contains with the from and to indices of the
neighbour links and their weights. `df2sn` converts the data retrieved
from reading the data from `df` back into a `spatial.neighbour` object.

## Author

Roger Bivand <Roger.Bivand@nhh.no>

## Note

Original idea due to Gidske Leknes Andersen, Department of Biology,
University of Bergen, Norway

## See also

[`sn2listw`](https://r-spatial.github.io/spdep/reference/listw2sn.md)

## Examples

``` r
columbus <- st_read(system.file("shapes/columbus.gpkg", package="spData")[1], quiet=TRUE)
col.gal.nb <- read.gal(system.file("weights/columbus.gal", package="spData")[1])
res <- listw2lines(nb2listw(col.gal.nb), st_geometry(columbus))
summary(res)
#>        i               j             i_ID               j_ID          
#>  Min.   : 1.00   Min.   : 1.00   Length:230         Length:230        
#>  1st Qu.:13.00   1st Qu.:13.00   Class :character   Class :character  
#>  Median :24.00   Median :24.00   Mode  :character   Mode  :character  
#>  Mean   :24.19   Mean   :24.19                                        
#>  3rd Qu.:35.00   3rd Qu.:35.00                                        
#>  Max.   :49.00   Max.   :49.00                                        
#>        wt               geometry  
#>  Min.   :0.1000   LINESTRING:230  
#>  1st Qu.:0.1429   epsg:NA   :  0  
#>  Median :0.1667                   
#>  Mean   :0.2130                   
#>  3rd Qu.:0.2500                   
#>  Max.   :0.5000                   
tf <- paste0(tempfile(), ".gpkg")
st_write(res, dsn=tf, driver="GPKG")
#> Writing layer `file280ec2a49b68e' to data source 
#>   `/tmp/Rtmp3LNfdt/file280ec2a49b68e.gpkg' using driver `GPKG'
#> Writing 230 features with 5 fields and geometry type Line String.
inMap <- st_read(tf)
#> Reading layer `file280ec2a49b68e' from data source 
#>   `/tmp/Rtmp3LNfdt/file280ec2a49b68e.gpkg' using driver `GPKG'
#> Simple feature collection with 230 features and 5 fields
#> Geometry type: LINESTRING
#> Dimension:     XY
#> Bounding box:  xmin: 6.165913 ymin: 11.04088 xmax: 10.96206 ymax: 14.43766
#> Projected CRS: Undefined Cartesian SRS with unknown unit
summary(inMap)
#>        i               j             i_ID               j_ID          
#>  Min.   : 1.00   Min.   : 1.00   Length:230         Length:230        
#>  1st Qu.:13.00   1st Qu.:13.00   Class :character   Class :character  
#>  Median :24.00   Median :24.00   Mode  :character   Mode  :character  
#>  Mean   :24.19   Mean   :24.19                                        
#>  3rd Qu.:35.00   3rd Qu.:35.00                                        
#>  Max.   :49.00   Max.   :49.00                                        
#>        wt                 geom    
#>  Min.   :0.1000   LINESTRING:230  
#>  1st Qu.:0.1429   epsg:NA   :  0  
#>  Median :0.1667                   
#>  Mean   :0.2130                   
#>  3rd Qu.:0.2500                   
#>  Max.   :0.5000                   
diffnb(sn2listw(df2sn(as.data.frame(inMap)))$neighbours, col.gal.nb)
#> Warning: style is M (missing); style should be set to a valid value
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
res1 <- listw2lines(nb2listw(col.gal.nb), as(columbus, "Spatial"))
summary(res1)
#> Object of class SpatialLinesDataFrame
#> Coordinates:
#>         min      max
#> x  6.221943 10.95359
#> y 11.010031 14.36908
#> Is projected: TRUE 
#> proj4string : [NA]
#> Data attributes:
#>        i               j             i_ID               j_ID          
#>  Min.   : 1.00   Min.   : 1.00   Length:230         Length:230        
#>  1st Qu.:13.00   1st Qu.:13.00   Class :character   Class :character  
#>  Median :24.00   Median :24.00   Mode  :character   Mode  :character  
#>  Mean   :24.19   Mean   :24.19                                        
#>  3rd Qu.:35.00   3rd Qu.:35.00                                        
#>  Max.   :49.00   Max.   :49.00                                        
#>        wt        
#>  Min.   :0.1000  
#>  1st Qu.:0.1429  
#>  Median :0.1667  
#>  Mean   :0.2130  
#>  3rd Qu.:0.2500  
#>  Max.   :0.5000  
```
