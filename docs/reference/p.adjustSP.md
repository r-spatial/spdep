# Adjust local association measures' p-values

Make an adjustment to local association measures' p-values based on the
number of neighbours (+1) of each region, rather than the total number
of regions.

## Usage

``` r
p.adjustSP(p, nb, method = "none")
```

## Arguments

- p:

  vector of p-values

- nb:

  a list of neighbours of class `nb`

- method:

  correction method as defined in
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html): "The adjustment
  methods include the Bonferroni correction ('"bonferroni"') in which
  the p-values are multiplied by the number of comparisons. Four less
  conservative corrections are also included by Holm (1979) ('holm'),
  Hochberg (1988) ('hochberg'), Hommel (1988) ('hommel') and Benjamini &
  Hochberg (1995) ('fdr'), respectively. A pass-through option ('none')
  is also included."

## Value

A vector of corrected p-values using only the number of neighbours + 1.

## Author

Danlin Yu and Roger Bivand <Roger.Bivand@nhh.no>

## See also

[`p.adjust`](https://rdrr.io/r/stats/p.adjust.html),
[`localG`](https://r-spatial.github.io/spdep/reference/localG.md),
[`localmoran`](https://r-spatial.github.io/spdep/reference/localmoran.md)

## Examples

``` r
data(afcon, package="spData")
oid <- order(afcon$id)
resG <- as.vector(localG(afcon$totcon, nb2listw(include.self(paper.nb))))
non <- format.pval(pnorm(2*(abs(resG)), lower.tail=FALSE), 2)
bon <- format.pval(p.adjustSP(pnorm(2*(abs(resG)), lower.tail=FALSE),
 paper.nb, "bonferroni"), 2)
tot <- format.pval(p.adjust(pnorm(2*(abs(resG)), lower.tail=FALSE),
 "bonferroni", n=length(resG)), 2)
data.frame(resG, non, bon, tot, row.names=afcon$name)[oid,]
#>                                 resG     non     bon     tot
#> THE GAMBIA               -0.98383592 0.02455 0.04911  1.0000
#> MALI                     -1.69893391 0.00034 0.00272  0.0143
#> SENEGAL                  -1.46321911 0.00171 0.00857  0.0720
#> BENIN                    -1.30139679 0.00462 0.02312  0.1942
#> MAURITANIA               -0.60496775 0.11315 0.56576  1.0000
#> NIGER                    -1.04877003 0.01797 0.14378  0.7549
#> IVORY COAST              -1.41712454 0.00230 0.01378  0.0965
#> GUINEA                   -1.44888005 0.00188 0.01128  0.0789
#> BURKINA FASO             -1.75085492 0.00023 0.00162  0.0097
#> LIBERIA                  -1.04053617 0.01871 0.07485  0.7860
#> SIERRA LEONE             -0.87032623 0.04087 0.12262  1.0000
#> GHANA                    -1.10269327 0.01371 0.05485  0.5760
#> TOGO                     -0.99053008 0.02379 0.09517  0.9993
#> CAMEROON                 -1.13328519 0.01171 0.07025  0.4917
#> NIGERIA                  -1.17261672 0.00951 0.04754  0.3993
#> GABON                    -0.78935857 0.05720 0.17160  1.0000
#> CENTRAL AFRICAN REPUBLIC  1.17349763 0.00946 0.05678  0.3974
#> CHAD                      0.46259185 0.17744 1.00000  1.0000
#> CONGO                    -0.20253005 0.34272 1.00000  1.0000
#> ZAIRE                     2.02270432 2.6e-05 0.00026  0.0011
#> ANGOLA                    1.23450728 0.00677 0.02710  0.2845
#> UGANDA                    3.33600851 1.3e-11 7.6e-11 5.3e-10
#> KENYA                     3.50301896 1.2e-12 7.4e-12 5.1e-11
#> TANZANIA                  1.09843592 0.01401 0.12613  0.5886
#> BURUNDI                   0.77417084 0.06077 0.24308  1.0000
#> RWANDA                    1.45720776 0.00178 0.00891  0.0748
#> SOMALIA                   1.18316273 0.00898 0.02695  0.3773
#> ETHIOPIA                  2.62720027 7.4e-08 3.0e-07 3.1e-06
#> ZAMBIA                    0.75273285 0.06610 0.59492  1.0000
#> ZIMBABWE                 -0.19956472 0.34490 1.00000  1.0000
#> MALAWI                    0.21195283 0.33582 1.00000  1.0000
#> MOZAMBIQUE               -0.28761679 0.28257 1.00000  1.0000
#> SOUTH AFRICA             -0.86814954 0.04126 0.33004  1.0000
#> LESOTHO                  -0.29841469 0.27531 0.55062  1.0000
#> BOTSWANA                  0.04090396 0.46740 1.00000  1.0000
#> SWAZILAND                -0.65938417 0.09362 0.28087  1.0000
#> MOROCCO                   0.02191606 0.48252 1.00000  1.0000
#> ALGERIA                  -0.36307938 0.23387 1.00000  1.0000
#> TUNISIA                   0.57910139 0.12339 0.37017  1.0000
#> LIBYA                     2.55272169 1.7e-07 1.2e-06 6.9e-06
#> SUDAN                     4.03925235 3.3e-16 3.0e-15 1.4e-14
#> EGYPT                     4.42133637 < 2e-16 < 2e-16 < 2e-16
```
