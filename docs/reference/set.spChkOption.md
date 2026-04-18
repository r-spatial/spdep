# Control checking of spatial object IDs

Provides support for checking the mutual integrity of spatial neighbour
weights and spatial data; similar mechanisms are used for passing global
verbose and zero.policy options, and for causing functions creating
neighbour objects to warn if there are multiple subgraphs.

## Usage

``` r
set.spChkOption(check)
get.spChkOption()
chkIDs(x, listw)
spNamedVec(var, data)
set.VerboseOption(check)
get.VerboseOption()
set.ZeroPolicyOption(check)
get.ZeroPolicyOption()
set.SubgraphOption(check)
get.SubgraphOption()
set.SubgraphCeiling(value)
get.SubgraphCeiling()
set.NoNeighbourOption(check)
get.NoNeighbourOption()
set.listw_is_CsparseMatrix_Option(check)
get.listw_is_CsparseMatrix_Option()
```

## Arguments

- check:

  a logical value, TRUE or FALSE

- value:

  an integer value, initialised as 100000L, the sum of the numbers of
  nodes and edges in the neighbour graph

- x:

  a vector the same length, or a two-dimensional array, or data frame
  with the same number of rows as the neighbours list in listw

- listw:

  a `listw` object or `nb` object inheriting from "nb"

- var:

  a character string or integer value for the column to be selected

- data:

  a two-dimensional array or data frame containing var

## Details

Analysis functions will have an spChk argument by default set to NULL,
and will call `get.spChkOption()` to get the global spatial option for
whether to check or not — this is initialised to FALSE, and consequently
should not break anything. It can be changed to TRUE using
`set.spChkOption(TRUE)`, or the spChk argument can be assigned in
analysis functions. `spNamedVec()` is provided to ensure that rownames
are passed on to single columns taken from two-dimensional arrays and
data frames.

## Value

`set.spChkOption()` returns the old logical value, `get.spChkOption()`
returns the current logical value, and `chkIDs()` returns a logical
value for the test lack of difference. `spNamedVec()` returns the
selected column with the names set to the row names of the object from
which it has been extracted.

## Author

Roger Bivand <Roger.Bivand@nhh.no>

## Note

The motivation for this mechanism is provided by the observation that
spatial objects on a map and their attribute data values need to be
linked uniquely, to avoid spurious results. The reordering between the
legacy Columbus data set used the earlier publications and that
available for download from the Spacestat website is just one example of
a common problem.

## Examples

``` r
data(oldcol)
rownames(COL.OLD)
#>  [1] "1001" "1002" "1003" "1004" "1005" "1006" "1007" "1008" "1009" "1010"
#> [11] "1011" "1012" "1013" "1014" "1015" "1016" "1017" "1018" "1019" "1020"
#> [21] "1021" "1022" "1023" "1024" "1025" "1026" "1027" "1028" "1029" "1030"
#> [31] "1031" "1032" "1033" "1034" "1035" "1036" "1037" "1038" "1039" "1040"
#> [41] "1041" "1042" "1043" "1044" "1045" "1046" "1047" "1048" "1049"
data(columbus, package="spData")
rownames(columbus)
#>  [1] "1005" "1001" "1006" "1002" "1007" "1008" "1004" "1003" "1018" "1010"
#> [11] "1038" "1037" "1039" "1040" "1009" "1036" "1011" "1042" "1041" "1017"
#> [21] "1043" "1019" "1012" "1035" "1032" "1020" "1021" "1031" "1033" "1034"
#> [31] "1045" "1013" "1022" "1044" "1023" "1046" "1030" "1024" "1047" "1016"
#> [41] "1014" "1049" "1029" "1025" "1028" "1048" "1015" "1027" "1026"
get.spChkOption()
#> [1] FALSE
oldChk <- set.spChkOption(TRUE)
get.spChkOption()
#> [1] TRUE
chkIDs(COL.OLD, nb2listw(COL.nb))
#> [1] TRUE
chkIDs(columbus, nb2listw(col.gal.nb))
#> [1] TRUE
chkIDs(columbus, nb2listw(COL.nb))
#> [1] FALSE
tmp <- try(moran.test(spNamedVec("CRIME", COL.OLD), nb2listw(COL.nb)))
tmp <- try(moran.test(spNamedVec("CRIME", columbus), nb2listw(col.gal.nb)))
tmp <- try(moran.test(spNamedVec("CRIME", columbus), nb2listw(COL.nb)))
#> Error in moran.test(spNamedVec("CRIME", columbus), nb2listw(COL.nb)) : 
#>   Check of data and weights ID integrity failed
set.spChkOption(FALSE)
get.spChkOption()
#> [1] FALSE
moran.test(spNamedVec("CRIME", columbus), nb2listw(COL.nb))
#> 
#>  Moran I test under randomisation
#> 
#> data:  spNamedVec("CRIME", columbus)  
#> weights: nb2listw(COL.nb)    
#> 
#> Moran I statistic standard deviate = 3.8402, p-value = 6.147e-05
#> alternative hypothesis: greater
#> sample estimates:
#> Moran I statistic       Expectation          Variance 
#>       0.341628707      -0.020833333       0.008908762 
#> 
tmp <- try(moran.test(spNamedVec("CRIME", columbus), nb2listw(COL.nb),
 spChk=TRUE), silent=TRUE)
set.spChkOption(oldChk)
get.spChkOption()
#> [1] FALSE
```
