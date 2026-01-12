# Measure distance from plot

Measure a distance between two points on a plot using `locator`; the
function checks `par("plt")` and `par("usr")` to try to ensure that the
aspect ratio y/x is 1, that is that the units of measurement in both x
and y are equivalent.

## Usage

``` r
airdist(ann=FALSE)
```

## Arguments

- ann:

  annotate the plot with line measured and distance

## Value

a list with members:

- dist:

  distance measured

- coords:

  coordinates between which distance is measured

## Author

Roger Bivand <Roger.Bivand@nhh.no>

## See also

[`locator`](https://rdrr.io/r/graphics/locator.html)
