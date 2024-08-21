% Copyright 2003-4 by Roger S. Bivand
\name{read.gwt2nb}
\alias{read.gwt2nb}
\alias{write.sn2gwt}
\alias{read.dat2listw}
\alias{read.swmdbf2listw}
\alias{read_swm_dbf}
\alias{write.sn2dat}
\alias{write.sn2DBF}
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
write.sn2DBF(sn, file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file, fn}{name of file with weights data}
  \item{region.id}{a character vector of region IDs - for ArcGIS SWM DBFs, the values must be character integers (only numbers not starting with zero)}
  \item{sn}{a \code{spatial.neighbour} object}
  \item{shpfile}{character string: if not given Shapefile name taken from GWT file for this dataset}
  \item{ind}{character string: region id indicator field name}
  \item{useInd}{default FALSE, if TRUE, write \code{region.id} attribute ID key tags to output file (use in OpenGeoDa will depend on the shapefile having the field named in the \code{ind} argument matching the exported tags)}
  \item{legacy}{default FALSE; if TRUE, header has single field with number of observations only}
  \item{style}{default NULL, missing, set to "M" and warning given; if not "M", passed to \code{\link{nb2listw}} to re-build the object}
  \item{zero.policy}{default NULL, use global option value; if FALSE stop with error for any empty neighbour sets, if TRUE permit the weights list to be formed with zero-length weights vectors}
}
\details{
Attempts to honour the region.id argument given when reading GWT and SWM/DBF files. If the region IDs given in \code{region.id=} do not match the origins or destinations in the GWT file, an error or warning will be thrown, which should be considered carefully. \code{read_swm_dbf} is a simplified interface to \code{read.swmdbf2listw} When no-neighbour observations are present, it is advised that \code{region.id=} be given in \code{read.swmdbf2listw}; while the function should read correctly if the minimum and maximum IDs are present as observations with neighbours without \code{region.id=} set, reading using \code{read.swmdbf2listw} will fail when the minimum or maximum ID observations have no neighbours and \code{region.id=} is not given.
}
\value{
  \code{read.gwt2nb} returns a neighbour "nb" object with the generalised weights stored as a list element called "dlist" of the "GeoDa" attribute; \code{read.swmdbf2listw} returns a "listw" object read from a DBF file exported from an ArcGIS SWM object.
}

\references{Luc Anselin (2003) \emph{GeoDa 0.9 User's Guide}, pp. 80--81, Spatial Analysis Laboratory, Department of Agricultural and Consumer Economics, University of Illinois, Urbana-Champaign, \url{http://geodacenter.github.io/docs/geoda093.pdf}; also material formerly at spatial-econometrics.com/data/contents.html}

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
if (require("foreign", quietly=TRUE)) {
nc_sf <- sf::st_read(system.file("gpkg/nc.gpkg", package="sf")[1])
nc_sf$UniqueID <- 1:nrow(nc_sf)
fn <- system.file("etc/misc/nc_contiguity_unique_id.dbf", package="spdep")[1]
nc1 <- read.swmdbf2listw(fn, style="B")
nc1
nc1a <- read.swmdbf2listw(fn, region.id=as.character(nc_sf$UniqueID),
 style="B")
all.equal(nc1, nc1a)
fn <- system.file("etc/misc/nc_contiguity_unique_id_islands.dbf",
 package="spdep")[1]
try(nc1i <- read.swmdbf2listw(fn, style="B"))
nc1i <- read.swmdbf2listw(fn, style="B", zero.policy=TRUE)
nc1ia <- read.swmdbf2listw(fn, region.id=as.character(nc_sf$UniqueID),
 style="B", zero.policy=TRUE)
nc1ia
all.equal(nc1i, nc1ia)
cal <- st_read(system.file("etc/shapes/california.gpkg", package="spdep")[1])
fn <- system.file("etc/misc/contiguity_myid.dbf", package="spdep")[1]
cal1 <- read.swmdbf2listw(fn, style="B")
cal1a <- read.swmdbf2listw(fn, region.id=as.character(cal$MYID), style="B")
all.equal(cal1, cal1a)
fn <- system.file("etc/misc/contiguity_unique_id.dbf", package="spdep")[1]
cal2 <- read.swmdbf2listw(fn, style="B")
cal2a <- read.swmdbf2listw(fn, region.id=as.character(cal$UniqueID), style="B")
all.equal(cal2, cal2a)
fn <- system.file("etc/misc/contiguity_unique_id_islands.dbf", package="spdep")[1]
try(cal3i <- read.swmdbf2listw(fn, style="B"))
cal3i <- read.swmdbf2listw(fn, style="B", zero.policy=TRUE)
cal3ia <- read.swmdbf2listw(fn, region.id=as.character(cal$UniqueID), style="B", zero.policy=TRUE)
all.equal(cal3i, cal3ia)
cal1a_1n_nb <- cal1a$neighbours
cal1a_1n_nb <- droplinks(cal1a_1n_nb, drop=c("158", "180", "215"), sym=TRUE)
cal1a_1n <- nb2listw(cal1a_1n_nb, style="B", zero.policy=TRUE)
cal1a_1n_sn <- listw2sn(cal1a_1n)
file <- tempfile(fileext=".dbf")
write.sn2DBF(cal1a_1n_sn, file)
cal1a_1n_rt <- read.swmdbf2listw(file, region.id=as.character(cal$MYID),
 style="B", zero.policy=TRUE)
all.equal(cal1a_1n$neighbours, cal1a_1n_rt$neighbours)
all.equal(cal1a_1n$weights, cal1a_1n_rt$weights, check.attributes=FALSE)
cal1_1n_rt <- read.swmdbf2listw(file, style="B", zero.policy=TRUE)
all(isTRUE(all.equal(cal1a_1n$neighbours, cal1_1n_rt$neighbours)))
all(isTRUE(all.equal(cal1a_1n$weights, cal1_1n_rt$weights)))
}
}
\keyword{spatial}
