# Prune a Minimun Spanning Tree

This function deletes a first edge and makes two subsets of edges. Each
subset is a Minimun Spanning Treee.

## Usage

``` r
prunemst(edges, only.nodes = TRUE)
```

## Arguments

- edges:

  A matrix with two colums with each row is one edge

- only.nodes:

  If `only.nodes=FALSE`, return a edges and nodes of each MST resulted.
  If `only.nodes=TRUE`, return a two sets of nodes. Defalt is TRUE

## Value

A list of length two. If `only.nodes=TRUE` each element is a vector of
nodes. If `only.nodes=FALSE` each element is a list with nodes and
edges.

## Author

Elias T. Krainski and Renato M. Assuncao

## See also

See Also as
[`mstree`](https://r-spatial.github.io/spdep/reference/mstree.md)

## Examples

``` r
e <- matrix(c(2,3, 1,2, 3,4, 4,5), ncol=2, byrow=TRUE)
e
#>      [,1] [,2]
#> [1,]    2    3
#> [2,]    1    2
#> [3,]    3    4
#> [4,]    4    5
prunemst(e)
#> $node1
#> [1] 2 1
#> 
#> $node2
#> [1] 3 4 5
#> 
prunemst(e, only.nodes=FALSE)
#> [[1]]
#> [[1]]$node
#> [1] 2 1
#> 
#> [[1]]$edge
#>      [,1] [,2]
#> [1,]    1    2
#> 
#> 
#> [[2]]
#> [[2]]$node
#> [1] 3 4 5
#> 
#> [[2]]$edge
#>      [,1] [,2]
#> [1,]    3    4
#> [2,]    4    5
#> 
#> 
```
