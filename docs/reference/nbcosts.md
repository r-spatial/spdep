# Compute cost of edges

The cost of each edge is the distance between it nodes. This function
compute this distance using a data.frame with observations vector in
each node.

## Usage

``` r
nbcost(data, id, id.neigh,  method = c("euclidean", "maximum", 
    "manhattan", "canberra", "binary", "minkowski", "mahalanobis"),
    p = 2, cov, inverted = FALSE)
nbcosts(nb, data,  method = c("euclidean", "maximum", 
    "manhattan", "canberra", "binary", "minkowski", "mahalanobis"),
    p = 2, cov, inverted = FALSE)
```

## Arguments

- nb:

  An object of `nb` class. See
  [`poly2nb`](https://r-spatial.github.io/spdep/reference/poly2nb.md)
  for details.

- data:

  A matrix with observations in the nodes.

- id:

  Node index to compute the cost

- id.neigh:

  Idex of neighbours nodes of node `id`

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

A object of `nbdist` class. See
[`nbdists`](https://r-spatial.github.io/spdep/reference/nbdists.md) for
details.

## Note

The neighbours must be a connected graph.

## Author

Elias T. Krainski and Renato M. Assuncao

## See also

See Also as
[`nbdists`](https://r-spatial.github.io/spdep/reference/nbdists.md),
[`nb2listw`](https://r-spatial.github.io/spdep/reference/nb2listw.md)
