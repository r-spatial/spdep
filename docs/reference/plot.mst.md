# Plot the Minimum Spanning Tree

This function plots a MST, the nodes are circles and the edges are
segments.

## Usage

``` r
# S3 method for class 'mst'
plot(x, coords, label.areas = NULL, 
   cex.circles = 1, cex.labels = 1, add=FALSE, ...)
```

## Arguments

- x:

  Object of `mst` class.

- coords:

  A two column matrix with the coordinates of nodes.

- label.areas:

  A vector with the labels of nodes

- cex.circles:

  The length of circles to plot.

- cex.labels:

  The length of nodes labels ploted.

- add:

  default FALSE, create new plot

- ...:

  Further arguments passed to plotting functions.

## Author

Elias T. Krainski and Renato M. Assuncao

## See also

See Also as
[`skater`](https://r-spatial.github.io/spdep/reference/skater.md) and
[`mstree`](https://r-spatial.github.io/spdep/reference/mstree.md)

## Examples

``` r
### see example in mstree function documentation
```
