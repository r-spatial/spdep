# Subset a spatial weights list

The function subsets a spatial weights list, retaining objects for which
the subset argument vector is TRUE. At present it will only subset
non-general weights lists (that is those created by `nb2listw` with
`glist=NULL`).

## Usage

``` r
# S3 method for class 'listw'
subset(x, subset, zero.policy = attr(x, "zero.policy"), ...)
```

## Arguments

- x:

  an object of class `listw`

- subset:

  logical expression

- zero.policy:

  default `attr(x, "zero.policy")` as set when `x` was created, if
  attribute not set, use global option value; if FALSE stop with error
  for any empty neighbour sets, if TRUE permit the weights list to be
  formed with zero-length weights vectors - passed through to `nb2listw`

- ...:

  generic function pass-through

## Value

The function returns an object of class `listw` with component `style`
the same as the input object, component `neighbours` a list of integer
vectors containing neighbour region number ids (compacted to run from
1:number of regions in subset), and component `weights` as the weights
computed for `neighbours` using `style`. If no-neighbour observations
are created by subsetting and `zero.policy` in the input weights object
was FALSE, it will be set to TRUE and a warning issued.

## Author

Roger Bivand <Roger.Bivand@nhh.no>

## See also

[`nb2listw`](https://r-spatial.github.io/spdep/reference/nb2listw.md),
[`subset.nb`](https://r-spatial.github.io/spdep/reference/subset.nb.md)

## Examples

``` r
col.gal.nb <- read.gal(system.file("weights/columbus.gal", package="spData")[1])
to.be.dropped <- c(31, 34, 36, 39, 42, 46)
pre <- nb2listw(col.gal.nb)
print(pre)
#> Characteristics of weights list object:
#> Neighbour list object:
#> Number of regions: 49 
#> Number of nonzero links: 230 
#> Percentage nonzero weights: 9.579342 
#> Average number of links: 4.693878 
#> 
#> Weights style: W 
#> Weights constants summary:
#>    n   nn S0       S1       S2
#> W 49 2401 49 23.48489 204.6687
post <- subset(pre, !(1:length(col.gal.nb) %in% to.be.dropped))
print(post)
#> Characteristics of weights list object:
#> Neighbour list object:
#> Number of regions: 43 
#> Number of nonzero links: 212 
#> Percentage nonzero weights: 11.46566 
#> Average number of links: 4.930233 
#> 
#> Weights style: W 
#> Weights constants summary:
#>    n   nn S0       S1       S2
#> W 43 1849 43 19.26584 178.4604
columbus <- st_read(system.file("shapes/columbus.gpkg", package="spData")[1], quiet=TRUE)
nb <- poly2nb(columbus)
lw <- nb2listw(nb, style="W")
attr(lw, "zero.policy")
#> [1] FALSE
(lwa <- subset(lw, 1:nrow(columbus) != c(21)))
#> Warning: subsetting caused increase in subgraph count
#> Characteristics of weights list object:
#> Neighbour list object:
#> Number of regions: 48 
#> Number of nonzero links: 230 
#> Percentage nonzero weights: 9.982639 
#> Average number of links: 4.791667 
#> 2 disjoint connected subgraphs
#> 
#> Weights style: W 
#> Weights constants summary:
#>    n   nn S0       S1       S2
#> W 48 2304 48 22.46811 199.4398
attr(lwa, "zero.policy")
#> [1] FALSE
(lwb <- subset(lw, !(1:nrow(columbus) %in% c(21, 36, 39))))
#> Warning: subsetting caused increase in subgraph count
#> Warning: subsetting created no-neighbour observations, zero.policy set TRUE
#> Characteristics of weights list object:
#> Neighbour list object:
#> Number of regions: 46 
#> Number of nonzero links: 216 
#> Percentage nonzero weights: 10.20794 
#> Average number of links: 4.695652 
#> 1 region with no links:
#> 46
#> 3 disjoint connected subgraphs
#> 
#> Weights style: W 
#> Weights constants summary:
#>    n   nn S0     S1       S2
#> W 45 2025 45 22.857 187.4843
attr(lwb, "zero.policy")
#> [1] TRUE
```
