# Columbus OH spatial analysis data set - old numbering

The `COL.OLD` data frame has 49 rows and 22 columns. The observations
are ordered and numbered as in the original analyses of the data set in
the SpaceStat documentation and in Anselin, L. 1988 Spatial
econometrics: methods and models, Dordrecht: Kluwer. Unit of analysis:
49 neighbourhoods in Columbus, OH, 1980 data. In addition the data set
includes `COL.nb`, the neighbours list as used in Anselin (1988).

## Usage

``` r
data(oldcol)
```

## Format

This data frame contains the following columns:

- AREA_PL:

  computed by ArcView (agrees with areas of polygons in the “columbus”
  data set

- PERIMETER:

  computed by ArcView

- COLUMBUS.:

  internal polygon ID (ignore)

- COLUMBUS.I:

  another internal polygon ID (ignore)

- POLYID:

  yet another polygon ID

- NEIG:

  neighborhood id value (1-49); conforms to id value used in Spatial
  Econometrics book.

- HOVAL:

  housing value (in \$1,000)

- INC:

  household income (in \$1,000)

- CRIME:

  residential burglaries and vehicle thefts per thousand households in
  the neighborhood

- OPEN:

  open space in neighborhood

- PLUMB:

  percentage housing units without plumbin

- DISCBD:

  distance to CBD

- X:

  x coordinate (in arbitrary digitizing units, not polygon coordinates)

- Y:

  y coordinate (in arbitrary digitizing units, not polygon coordinates)

- AREA_SS:

  neighborhood area (computed by SpaceStat)

- NSA:

  north-south dummy (North=1)

- NSB:

  north-south dummy (North=1)

- EW:

  east-west dummy (East=1)

- CP:

  core-periphery dummy (Core=1)

- THOUS:

  constant=1,000

- NEIGNO:

  NEIG+1,000, alternative neighborhood id value

- PERIM:

  polygon perimeter (computed by SpaceStat)

## Details

The row names of `COL.OLD` and the `region.id` attribute of `COL.nb` are
set to `columbus$NEIGNO`.

## Source

Anselin, Luc. 1988. Spatial econometrics: methods and models. Dordrecht:
Kluwer Academic, Table 12.1 p. 189.

## Note

All source data files prepared by Luc Anselin, Spatial Analysis
Laboratory, Department of Agricultural and Consumer Economics,
University of Illinois, Urbana-Champaign.
