# Compute the sum of dissimilarity

This function computes the sum of dissimilarity between each observation
and the mean (scalar of vector) of the observations.

## Usage

``` r
ssw(data, id, method = c("euclidean", "maximum", 
    "manhattan", "canberra", "binary", "minkowski",
    "mahalanobis"), p = 2, cov, inverted = FALSE)
```

## Arguments

- data:

  A matrix with observations in the nodes.

- id:

  Node index to compute the cost

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

A numeric, the sum of dissimilarity between the observations `id` of
`data` and the mean (scalar of vector) of this observations.

## Author

Elias T. Krainski and Renato M. Assuncao

## See also

See Also as
[`nbcost`](https://r-spatial.github.io/spdep/reference/nbcosts.md)

## Examples

``` r
data(USArrests)
n <- nrow(USArrests)
ssw(USArrests, 1:n)
#> [1] 3701.394
ssw(USArrests, 1:(n/2))
#> [1] 1910.214
ssw(USArrests, (n/2+1):n)
#> [1] 1625.882
ssw(USArrests, 1:(n/2)) + ssw(USArrests, (n/2+1):n)
#> [1] 3536.096
```
