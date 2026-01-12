# Cardinalities for neighbours lists

The function tallies the numbers of neighbours of regions in the
neighbours list.

## Usage

``` r
card(nb)
```

## Arguments

- nb:

  a neighbours list object of class `nb`

## Value

An integer vector of the numbers of neighbours of regions in the
neighbours list.

## Details

“nb” objects are stored as lists of integer vectors, where the vectors
contain either the indices in the range `1:n` for `n` as `length(nb)` of
the neighbours of region `i`, or `as.integer(0)` to signal no
neighbours. The function `card(nb)` is used to extract the numbers of
neighbours from the “nb” object.

## References

Bivand R, Pebesma EJ, Gomez-Rubio V, (2008) *Applied Spatial Data
Analysis with R*, Springer, New York, pp. 239-251; Bivand R, Portnov B,
(2004) Exploring spatial data analysis techniques using R: the case of
observations with no neighbours. In: Anselin L, Florax R, Rey S, (eds.),
*Advances in Spatial Econometrics, Methodology, Tools and Applications*.
Berlin: Springer-Verlag, pp. 121-142,
[doi:10.1007/978-3-662-05617-2_6](https://doi.org/10.1007/978-3-662-05617-2_6)
.

## Author

Roger Bivand <Roger.Bivand@nhh.no>

## See also

[`summary.nb`](https://r-spatial.github.io/spdep/reference/summary.nb.md)

## Examples

``` r
col.gal.nb <- read.gal(system.file("weights/columbus.gal", package="spData")[1])
table(card(col.gal.nb))
#> 
#>  2  3  4  5  6  7  8  9 10 
#>  7  7 13  4  9  6  1  1  1 
```
