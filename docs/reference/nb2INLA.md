# Output spatial neighbours for INLA

Output spatial neighbours for INLA

## Usage

``` r
nb2INLA(file, nb)
```

## Arguments

- file:

  file where adjacency matrix will be stored

- nb:

  an object of class `nb`

## Value

Nothing is returned but a file will be created with the representation
of the adjacency matrix as required by INLA for its spatial models.

## References

http://www.r-inla.org

## Author

Virgilio Gomez-Rubio

## Examples

``` r
col.gal.nb <- read.gal(system.file("weights/columbus.gal", package="spData")[1])
td <- tempdir()
x <- nb2INLA(paste(td, "columbus-INLA.adj", sep="/"), col.gal.nb)
readLines(paste(td, "columbus-INLA.adj", sep="/"), n=10)
#>  [1] "49"                        "1 2 2 3"                  
#>  [3] "2 3 1 3 4"                 "3 4 1 2 4 5"              
#>  [5] "4 4 2 3 5 8"               "5 7 3 4 6 8 9 11 15"      
#>  [7] "6 2 5 9"                   "7 4 8 12 13 14"           
#>  [9] "8 6 4 5 7 11 12 13"        "9 8 5 6 10 15 20 22 25 26"
```
