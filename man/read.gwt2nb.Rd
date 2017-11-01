% Copyright 2003-4 by Roger S. Bivand
\name{read.gwt2nb}
\alias{read.gwt2nb}
\alias{write.sn2gwt}
\alias{read.dat2listw}
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
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{name of file with weights data}
  \item{region.id}{region IDs}
  \item{sn}{a \code{spatial.neighbour} object}
  \item{shpfile}{character string: if not given Shapefile name taken from GWT file for this dataset}
  \item{ind}{character string: region id indicator field name}
  \item{useInd}{default FALSE, if TRUE, write \code{region.id} attribute ID key tags to output file (use in OpenGeoDa will depend on the shapefile having the field named in the \code{ind} argument matching the exported tags)}
  \item{legacy}{default FALSE; if TRUE, header has single field with number of observations only}
}
\details{
 Now attempts to honour the region.id argument given when reading GWT files. 
}
\value{
  \code{read.gwt2nb} returns a neighbour "nb" object with the generalised weights stored as a list element called "dlist" of the "GeoDa" attribute.
}

\references{Luc Anselin (2003) \emph{GeoDa 0.9 User's Guide}, pp. 80--81, Spatial Analysis Laboratory, Department of Agricultural and Consumer Economics, University of Illinois, Urbana-Champaign, \url{https://s3.amazonaws.com/geoda/software/docs/geoda093.pdf}; also \url{http://spatial-econometrics.com/data/contents.html}}

\author{Roger Bivand \email{Roger.Bivand@nhh.no}}
 
\seealso{\code{\link{read.gal}}}

\examples{
data(baltimore, package="spData")
STATION <- baltimore$STATION
gwt1 <- read.gwt2nb(system.file("etc/weights/baltk4.GWT", package="spdep")[1],
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
}
\keyword{spatial}
