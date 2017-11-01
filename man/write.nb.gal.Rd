% Copyright 2001-3 by Roger S. Bivand
\name{write.nb.gal}
\alias{write.nb.gal}

\title{Write a neighbours list as a GAL lattice file}
\description{
  Write a neighbours list as a GAL lattice file, may also use newer GeoDa header format
}
\usage{
write.nb.gal(nb, file, oldstyle=TRUE, shpfile=NULL, ind=NULL)
}

\arguments{
  \item{nb}{an object of class \code{nb} with a list of integer vectors containing neighbour region number ids.}
  \item{file}{name of file with GAL lattice data}
  \item{oldstyle}{if TRUE, first line of file contains only number of spatial units, if FALSE, uses newer GeoDa style}
  \item{shpfile}{Shapefile name taken from GAL file for this dataset}
  \item{ind}{region id indicator variable name}
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{read.gal}}}

\examples{
if (require(rgdal, quietly=TRUE)) {
example(columbus, package="spData")
GALfile <- tempfile("GAL")
write.nb.gal(col.gal.nb, GALfile)
col.queen <- read.gal(GALfile)
summary(diffnb(col.queen, col.gal.nb))
}}
\keyword{spatial}
