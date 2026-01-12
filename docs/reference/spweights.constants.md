# Provides constants for spatial weights matrices

The function calculates the constants needed for tests of spatial
autocorrelation for general weights matrices represented as `listw`
objects. Note: from spdep 0.3-32, the values of S1 and S2 are returned
correctly for both underlying symmetric and asymmetric neighbour lists,
before 0.3-32, S1 and S2 were wrong for listw objects based on
asymmetric neighbour lists, such as k-nearest neighbours (thanks to Luc
Anselin for finding the bug).

## Usage

``` r
spweights.constants(listw, zero.policy=attr(listw, "zero.policy"), adjust.n=TRUE)
Szero(listw)
```

## Arguments

- listw:

  a `listw` object from for example `nb2listw`

- zero.policy:

  default `attr(listw, "zero.policy")` as set when `listw` was created,
  if attribute not set, use global option value; if TRUE ignore zones
  without neighbours, if FALSE fail when encountered

- adjust.n:

  default TRUE, if FALSE the number of observations is not adjusted for
  no-neighbour observations, if TRUE, the number of observations is
  adjusted

## Value

- n:

  number of zones

- n1:

  n - 1

- n2:

  n - 2

- n3:

  n - 3

- nn:

  n \* n

- S0:

  global sum of weights

- S1:

  S1 sum of weights

- S2:

  S2 sum of weights

## References

Haining, R. 1990 Spatial data analysis in the social and environmental
sciences, Cambridge University Press, p. 233; Cliff, A. D., Ord, J. K.
1981 Spatial processes, Pion, p. 19, 21.

## Author

Roger Bivand <Roger.Bivand@nhh.no>

## See also

[`nb2listw`](https://r-spatial.github.io/spdep/reference/nb2listw.md)

## Examples

``` r
data(oldcol)
B <- spweights.constants(nb2listw(COL.nb, style="B"))
W <- spweights.constants(nb2listw(COL.nb, style="W"))
C <- spweights.constants(nb2listw(COL.nb, style="C"))
S <- spweights.constants(nb2listw(COL.nb, style="S"))
U <- spweights.constants(nb2listw(COL.nb, style="U"))
print(data.frame(rbind(unlist(B), unlist(W), unlist(C), unlist(S), unlist(U)),
  row.names=c("B", "W", "C", "S", "U")))
#>    n n1 n2 n3   nn  S0           S1           S2
#> B 49 48 47 46 2401 232 464.00000000 5.136000e+03
#> W 49 48 47 46 2401  49  23.29434146 2.048729e+02
#> C 49 48 47 46 2401  49  20.69827586 2.291085e+02
#> S 49 48 47 46 2401  49  21.25561347 2.134568e+02
#> U 49 48 47 46 2401   1   0.00862069 9.542212e-02
```
