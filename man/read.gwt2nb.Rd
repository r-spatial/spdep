% Copyright 2003-4 by Roger S. Bivand
\name{read.gwt2nb}
\alias{read.gwt2nb}
\alias{write.sn2gwt}
\alias{read.dat2listw}
\alias{read.swmdbf2listw}
\alias{read_swm_dbf}
\alias{write.sn2dat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Read and write spatial neighbour files}
\description{
  The "gwt" functions read and write GeoDa GWT files (the example file baltk4.GWT was downloaded from the site given in the reference), and the "dat" functions read and write Matlab sparse matrix files as used by James LeSage's Spatial Econometrics Toolbox (the example file wmat.dat was downloaded from the site given in the reference). The body of the files after any headers should have three columns separated by white space, and the third column must be numeric in the locale of the reading platform (correct decimal separator).
}
\usage{
read.gwt2nb(file, region.id=NULL)
write.sn2gwt(sn, file, shpfile=NULL, ind=NULL, useInd=FALSE, legacy=FALSE)
read.dat2listw(file)
write.sn2dat(sn, file)
read.swmdbf2listw(fn, region.id=NULL, style=NULL, zero.policy=NULL)
read_swm_dbf(fn)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file, fn}{name of file with weights data}
  \item{region.id}{region IDs}
  \item{sn}{a \code{spatial.neighbour} object}
  \item{shpfile}{character string: if not given Shapefile name taken from GWT file for this dataset}
  \item{ind}{character string: region id indicator field name}
  \item{useInd}{default FALSE, if TRUE, write \code{region.id} attribute ID key tags to output file (use in OpenGeoDa will depend on the shapefile having the field named in the \code{ind} argument matching the exported tags)}
  \item{legacy}{default FALSE; if TRUE, header has single field with number of observations only}
  \item{style}{default NULL, missing, set to "M" and warning given; if not "M", passed to \code{\link{nb2listw}} to re-build the object}
  \item{zero.policy}{default NULL, use global option value; if FALSE stop with error for any empty neighbour sets, if TRUE permit the weights list to be formed with zero-length weights vectors}
}
\details{
Attempts to honour the region.id argument given when reading GWT and SWM/DBF files. If the region IDs given in \code{region.id=} do not match the origins or destinations in the GWT file, an error or warning will be thrown, which should be considered carefully. \code{read_swm_dbf} is a simplified interface to \code{read.swmdbf2listw}.
}
\value{
  \code{read.gwt2nb} returns a neighbour "nb" object with the generalised weights stored as a list element called "dlist" of the "GeoDa" attribute; \code{read.swmdbf2listw} returns a "listw" object read from a DBF file exported from an ArcGIS SWM object.
}

\references{Luc Anselin (2003) \emph{GeoDa 0.9 User's Guide}, pp. 80--81, Spatial Analysis Laboratory, Department of Agricultural and Consumer Economics, University of Illinois, Urbana-Champaign, \url{http://geodacenter.github.io/docs/geoda093.pdf}; also \url{http://spatial-econometrics.com/data/contents.html}}

\author{Roger Bivand \email{Roger.Bivand@nhh.no}}
 
\seealso{\code{\link{read.gal}}}

\examples{
data(baltimore, package="spData")
STATION <- baltimore$STATION
gwt1 <- read.gwt2nb(system.file("weights/baltk4.GWT", package="spData")[1],
 STATION)
cat(paste("Neighbours list symmetry;", is.symmetric.nb(gwt1, FALSE, TRUE),
 "\n"))
listw1 <- nb2listw(gwt1, style="B", glist=attr(gwt1, "GeoDa")$dist)
tmpGWT <- tempfile()
write.sn2gwt(listw2sn(listw1), tmpGWT)
gwt2 <- read.gwt2nb(tmpGWT, STATION)
cat(paste("Neighbours list symmetry;", is.symmetric.nb(gwt2, FALSE, TRUE),
 "\n"))
diffnb(gwt1, gwt2)
data(oldcol)
tmpMAT <- tempfile()
COL.W <- nb2listw(COL.nb)
write.sn2dat(listw2sn(COL.W), tmpMAT)
listwmat1 <- read.dat2listw(tmpMAT)
diffnb(listwmat1$neighbours, COL.nb, verbose=TRUE)
listwmat2 <- read.dat2listw(system.file("etc/weights/wmat.dat", 
 package="spdep")[1])
diffnb(listwmat1$neighbours, listwmat2$neighbours, verbose=TRUE)
nc_sf <- sf::st_read(system.file("gpkg/nc.gpkg", package="sf")[1])
nc_sf$UniqueID <- 1:nrow(nc_sf)
fn <- system.file("etc/misc/nc_contiguity_unique_id.dbf", package="spdep")[1]
nc1 <- read.swmdbf2listw(fn, style="B")
nc1a <- read.swmdbf2listw(fn, region.id=as.character(nc_sf$UniqueID),
 style="B")
all.equal(nc1, nc1a)
fn <- system.file("etc/misc/nc_contiguity_unique_id_islands.dbf",
 package="spdep")[1]
try(nc1i <- read.swmdbf2listw(fn, style="B"))
nc1i <- read.swmdbf2listw(fn, style="B", zero.policy=TRUE)
nc1ia <- read.swmdbf2listw(fn, region.id=as.character(nc_sf$UniqueID),
 style="B", zero.policy=TRUE)
all.equal(nc1i, nc1ia)
}
\keyword{spatial}
