% Copyright 2001 by Roger S. Bivand
\name{oldcol}
\alias{oldcol}
\alias{COL.OLD}
\alias{COL.nb}
%\alias{coords.OLD}
%\alias{polys.OLD}
\docType{data}
\title{Columbus OH spatial analysis data set - old numbering}
\usage{data(oldcol)}
\description{
The \code{COL.OLD} data frame has 49 rows and 22 columns. The observations are
ordered and numbered as in the original analyses of the data set in the
SpaceStat documentation and in Anselin, L. 1988 Spatial econometrics: methods and models, Dordrecht: Kluwer.
Unit of analysis: 49 neighbourhoods in Columbus, OH, 1980 data.
In addition the data set includes \code{COL.nb}, the neighbours list
as used in Anselin (1988).
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{AREA_PL}{computed by ArcView (agrees with areas of polygons in the \dQuote{columbus} data set}
    \item{PERIMETER}{computed by ArcView}
    \item{COLUMBUS.}{internal polygon ID (ignore)}
    \item{COLUMBUS.I}{another internal polygon ID (ignore)}
    \item{POLYID}{yet another polygon ID}
    \item{NEIG}{neighborhood id value (1-49);
             conforms to id value used in Spatial Econometrics book.}
    \item{HOVAL}{housing value (in \$1,000)}
    \item{INC}{household income (in \$1,000)}
    \item{CRIME}{residential burglaries and vehicle thefts per thousand
        households in the neighborhood}
    \item{OPEN}{open space in neighborhood}
    \item{PLUMB}{percentage housing units without plumbin}
    \item{DISCBD}{distance to CBD}
    \item{X}{x coordinate (in arbitrary digitizing units, not polygon coordinates)}
    \item{Y}{y coordinate (in arbitrary digitizing units, not polygon coordinates)}
    \item{AREA_SS}{neighborhood area (computed by SpaceStat)}
    \item{NSA}{north-south dummy (North=1)}
    \item{NSB}{north-south dummy (North=1)}
    \item{EW}{east-west dummy (East=1)}
    \item{CP}{core-periphery dummy (Core=1)}
    \item{THOUS}{constant=1,000}
    \item{NEIGNO}{NEIG+1,000, alternative neighborhood id value}
    \item{PERIM}{polygon perimeter (computed by SpaceStat)}
  }
}
\details{
The row names of \code{COL.OLD} and the \code{region.id} attribute of
\code{COL.nb} are set to \code{columbus$NEIGNO}.
}
\source{
Anselin, Luc.  1988.  Spatial econometrics: methods and
models.  Dordrecht: Kluwer Academic, Table 12.1 p. 189.
}
\note{
All source data files prepared by
Luc Anselin, Spatial Analysis Laboratory, Department of Agricultural and Consumer Economics, University of Illinois, Urbana-Champaign, \url{https://spatial.uchicago.edu/sample-data}.
}
\keyword{datasets}
