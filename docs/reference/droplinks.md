# Drop and add links in a neighbours list

`droplinks` drops links to and from or just to a region from a
neighbours list. The example corresponds to Fingleton's Table 1, (1999)
p. 6, for lattices 5 to 19. `addlinks1` adds links from a single region
to specified regions.

## Usage

``` r
droplinks(nb, drop, sym=TRUE)
addlinks1(nb, from, to, sym=TRUE)
```

## Arguments

- nb:

  a neighbours list object of class `nb`

- drop:

  either a logical vector the length of `nb`, or a character vector of
  named regions corresponding to `nb`'s region.id attribute, or an
  integer vector of region numbers

- sym:

  TRUE for removal of both "row" and "column" links, FALSE for only
  "row" links; when adding links, inserts links to the from region from
  the to regions

- from:

  single from region for adding links, either a character vector of
  length 1 of the named from region corresponding to `nb`'s region.id
  attribute, or an integer vector of length 1 holding a region number

- to:

  to regions, either a character vector of named from regions
  corresponding to `nb`'s region.id attribute, or an integer vector of
  region numbers

## Value

The function returns an object of class `nb` with a list of integer
vectors containing neighbour region number ids.

## References

B. Fingleton (1999) Spurious spatial regression: some Monte Carlo
results with a spatial unit root and spatial cointegration, Journal of
Regional Science 39, pp. 1–19.

## Author

Roger Bivand <Roger.Bivand@nhh.no>

## See also

[`is.symmetric.nb`](https://r-spatial.github.io/spdep/reference/testnb.md)

## Examples

``` r
# \donttest{
rho <- c(0.2, 0.5, 0.95, 0.999, 1.0)
ns <- c(5, 7, 9, 11, 13, 15, 17, 19)
mns <- matrix(0, nrow=length(ns), ncol=length(rho))
rownames(mns) <- ns
colnames(mns) <- rho
mxs <- matrix(0, nrow=length(ns), ncol=length(rho))
rownames(mxs) <- ns
colnames(mxs) <- rho
for (i in 1:length(ns)) {
  nblist <- cell2nb(ns[i], ns[i])
  nbdropped <- droplinks(nblist, ((ns[i]*ns[i])+1)/2, sym=FALSE)
  listw <- nb2listw(nbdropped, style="W", zero.policy=TRUE)
  wmat <- listw2mat(listw)
  for (j in 1:length(rho)) {
    mat <- diag(ns[i]*ns[i]) - rho[j] * wmat
    res <- diag(solve(t(mat) %*% mat))
    mns[i,j] <- mean(res)
    mxs[i,j] <- max(res)
  }
}
#> Warning: some observations have no neighbours
#> Warning: some observations have no neighbours
#> Warning: some observations have no neighbours
#> Warning: some observations have no neighbours
#> Warning: some observations have no neighbours
#> Warning: some observations have no neighbours
#> Warning: some observations have no neighbours
#> Warning: some observations have no neighbours
print(mns)
#>         0.2      0.5      0.95     0.999          1
#> 5  1.038271 1.312627  9.486051  30.81487   32.04915
#> 7  1.036443 1.295621 10.899580  83.25437   92.09812
#> 9  1.035356 1.285145 10.798611 160.90951  195.02166
#> 11 1.034639 1.278279 10.383083 254.83998  347.71145
#> 13 1.034132 1.273442  9.968389 353.66366  555.88699
#> 15 1.033753 1.269852  9.619387 447.19245  824.46560
#> 17 1.033460 1.267082  9.337167 528.49015 1157.77630
#> 19 1.033227 1.264879  9.109487 594.23907 1559.69614
print(mxs)
#>         0.2      0.5     0.95     0.999          1
#> 5  1.048834 1.401934 12.00215  39.22742   40.79967
#> 7  1.048834 1.402174 14.66823 106.90031  118.01556
#> 9  1.048834 1.402176 15.49606 207.28928  249.74893
#> 11 1.048834 1.402176 15.75744 329.22973  443.97194
#> 13 1.048834 1.402176 15.83957 458.75739  707.14827
#> 15 1.048834 1.402176 15.86474 583.50722 1044.75562
#> 17 1.048834 1.402176 15.87225 695.10288 1461.57017
#> 19 1.048834 1.402176 15.87445 789.50575 1961.84025

# }
```
