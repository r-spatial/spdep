% Copyright 2001-3 by Roger S. Bivand
\name{read.gal}
\alias{read.gal}
\alias{read.geoda}
\title{Read a GAL lattice file into a neighbours list}
\description{
The function \code{read.gal()} reads a GAL lattice file into a neighbours list for spatial analysis. It will read old and new style (GeoDa) GAL files. The function \code{read.geoda} is a helper file for reading comma separated value data files, calling \code{read.csv()}.
}
\usage{
read.gal(file, region.id=NULL, override.id=FALSE)
read.geoda(file, row.names=NULL, skip=0)
}

\arguments{
  \item{file}{name of file with GAL lattice data}
  \item{region.id}{region IDs in specified order to coerse neighbours list order and numbering to that of the region.id}
  \item{override.id}{override any given (or NULL) region.id, collecting region.id numbering and order from the GAL file.}
  \item{row.names}{as in row.names in \code{read.csv()}, typically a character string naming the column of the file to be used}
  \item{skip}{skip number of lines, as in \code{read.csv()}}
}
\details{
Luc Anselin (2003): Spatial Analysis Laboratory, Department of Agricultural and Consumer Economics, University of Illinois, Urbana-Champaign, \url{http://www.csiss.org/gispopsci/workshops/2011/PSU/readings/W15_Anselin2007.pdf}; Luc Anselin (2003) \emph{GeoDa 0.9 User's Guide}, pp. 80--81, Spatial Analysis Laboratory, Department of Agricultural and Consumer Economics, University of Illinois, Urbana-Champaign, \url{https://s3.amazonaws.com/geoda/software/docs/geoda093.pdf}; GAL - Geographical Algorithms Library, University of Newcastle}

\value{
The function \code{read.gal()} returns an object of class \code{nb} with a list of integer vectors containing neighbour region number ids. The function \code{read.geoda} returns a data frame, and issues a warning if the returned object has only one column.
}


\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\note{Example data originally downloaded from now dead link: http://sal.agecon.uiuc.edu/weights/zips/us48.zip}

\seealso{\code{\link{summary.nb}}}

\examples{
us48.fipsno <- read.geoda(system.file("etc/weights/us48.txt",
 package="spdep")[1])
us48.q <- read.gal(system.file("etc/weights/us48_q.GAL", package="spdep")[1],
 us48.fipsno$Fipsno)
us48.r <- read.gal(system.file("etc/weights/us48_rk.GAL", package="spdep")[1],
 us48.fipsno$Fipsno)
data(state)
if (as.numeric(paste(version$major, version$minor, sep="")) < 19) {
 m50.48 <- match(us48.fipsno$"State.name", state.name)
} else {
 m50.48 <- match(us48.fipsno$"State_name", state.name)
}
plot(us48.q, as.matrix(as.data.frame(state.center))[m50.48,])
plot(diffnb(us48.r, us48.q),
 as.matrix(as.data.frame(state.center))[m50.48,], add=TRUE, col="red")
title(main="Differences between rook and queen criteria imported neighbours lists")
}
\keyword{spatial}
