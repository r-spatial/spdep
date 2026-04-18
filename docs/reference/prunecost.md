# Compute cost of prune each edge

If any edge are dropped, the MST are pruned. This generate a two
subgraphs. So, it makes a tree graphs and tree dissimilarity values are
computed, one for each graph. The dissimilarity is the sum over sqared
differences between the observactions in the nodes and mean vector of
observations in the graph. The dissimilarity of original graph and the
sum of dissimilarity of subgraphs are returned.

## Usage

``` r
prunecost(edges, data, method = c("euclidean", "maximum", "manhattan", 
    "canberra", "binary", "minkowski", "mahalanobis"), 
    p = 2, cov, inverted = FALSE)
```

## Arguments

- edges:

  A matrix with 2 colums with each row is one edge

- data:

  A data.frame with observations in the nodes.

- method:

  Character or function to declare distance method. If `method` is
  character, method must be "mahalanobis" or "euclidean", "maximum",
  "manhattan", "canberra", "binary" or "minkowisk". If `method` is one
  of "euclidean", "maximum", "manhattan", "canberra", "binary" or
  "minkowisk", see [`dist`](https://rdrr.io/r/stats/dist.html) for
  details, because this function as used to compute the distance. If
  `method="mahalanobis"`, the mahalanobis distance is computed between
  neighbour areas. If `method` is a `function`, this function is used to
  compute the distance.

- p:

  The power of the Minkowski distance.

- cov:

  The covariance matrix used to compute the mahalanobis distance.

- inverted:

  logical. If 'TRUE', 'cov' is supposed to contain the inverse of the
  covariance matrix.

## Value

A vector with the differences between the dissimilarity of all nodes and
the dissimilarity sum of all subgraphs obtained by pruning one edge each
time.

## Author

Elias T. Krainski and Renato M. Assuncao

## See also

See Also as
[`prunemst`](https://r-spatial.github.io/spdep/reference/prunemst.md)

## Examples

``` r
d <- data.frame(a=-2:2, b=runif(5))
e <- matrix(c(1,2, 2,3, 3,4, 4,5), ncol=2, byrow=TRUE)

sum(sweep(d, 2, colMeans(d))^2)
#> [1] 10.22026

prunecost(e, d)
#> [1] 2.161961 2.976434 3.197389 2.235287
```
