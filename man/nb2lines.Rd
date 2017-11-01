\name{nb2lines}
\alias{nb2lines}
\alias{listw2lines}
\alias{df2sn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Use arc-type shapefiles for import and export of weights}
\description{
  Use arc-type shapefiles for import and export of weights, storing spatial entity coordinates in the arcs, and the entity indices in the data frame.
}
\usage{
nb2lines(nb, wts, coords, proj4string=CRS(as.character(NA)))
listw2lines(listw, coords, proj4string=CRS(as.character(NA)))
df2sn(df, i="i", i_ID="i_ID", j="j", wt="wt")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nb}{a neighbour object of class \code{nb}}
  \item{wts}{list of general weights corresponding to neighbours}
  \item{coords}{matrix of region point coordinates}
  \item{proj4string}{Object of class CRS; holding a valid proj4 string}
  \item{listw}{a \code{listw} object of spatial weights}
  \item{df}{a data frame read from a shapefile, derived from the output of \code{nb2lines}}
  \item{i}{character name of column in df with from entity index}
  \item{i_ID}{character name of column in df with from entity region ID}
  \item{j}{character name of column in df with to entity index}
  \item{wt}{character name of column in df with weights}
}
\details{
 The maptools package function \code{writeSpatialShape} is used to transport out the list of lines made by \code{nb2lines} or \code{listw2lines}, which is a simple wrapper function. The neighbour and weights objects may be retrieved by converting the specified columns of the data slot of the SpatialLinesDataFrame object into a spatial.neighbour object, which is then converted into a weights list object.
}
\value{
  \code{nb2lines} and \code{listw2lines} return a SpatialLinesDataFrame object; its data slot contains a data frame with the from and to indices of the neighbour links and their weights. \code{df2sn} converts the data retrieved from reading the data from \code{df} back into a \code{spatial.neighbour} object.
}

\author{Roger Bivand \email{Roger.Bivand@nhh.no}}
\note{Original idea due to Gidske Leknes Andersen, Department of Biology, University of Bergen, Norway}

 
\seealso{\code{\link{sn2listw}}, \code{\link[maptools]{readShapeLines}}}

\examples{
require(maptools)
example(columbus)
coords <- coordinates(columbus)
res <- listw2lines(nb2listw(col.gal.nb), coords)
summary(res)
fn <- paste(tempdir(), "nbshape", sep="/")
writeLinesShape(res, fn=fn)
inMap <- readShapeLines(fn)
summary(inMap)
diffnb(sn2listw(df2sn(as(inMap, "data.frame")))$neighbours, col.gal.nb)
}
\keyword{spatial}

