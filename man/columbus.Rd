% Copyright 2017 by Roger S. Bivand
\name{columbus}
\alias{columbus}
\alias{col.gal.nb}
\alias{coords}
\alias{bbs}
\alias{polys}

\docType{data}
\title{Columbus OH spatial analysis data set}
\usage{data(columbus)}
\description{
The data set is now part of the spData package
}
\examples{
columbus <- st_read(system.file("shapes/columbus.gpkg", package="spData")[1], quiet=TRUE)
col.gal.nb <- read.gal(system.file("weights/columbus.gal", package="spData")[1])
}


\keyword{datasets}
