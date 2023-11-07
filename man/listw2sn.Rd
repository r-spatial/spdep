% Copyright 2001-3 by Roger S. Bivand
\name{listw2sn}
\alias{listw2sn}
\alias{sn2listw}
%\alias{as.spam.listw}
\title{Spatial neighbour sparse representation}
\description{
The function makes a \code{"spatial neighbour"} object representation
(similar to the S-PLUS spatial statististics module representation of a
\code{"listw"} spatial weights object. \code{sn2listw()} is the inverse function to \code{listw2sn()}, creating a \code{"listw"} object from a \code{"spatial neighbour"} object. }
\usage{
listw2sn(listw)
sn2listw(sn, style = NULL, zero.policy = NULL, from_mat2listw=FALSE)
}
\arguments{
  \item{listw}{a \code{listw} object from for example \code{nb2listw}}
  \item{sn}{a \code{spatial.neighbour} object}
  \item{style}{default NULL, missing, set to "M" and warning given; if not "M", passed to \code{\link{nb2listw}} to re-build the object}
  \item{zero.policy}{default NULL, use global option value; if FALSE stop with error for any empty neighbour sets, if TRUE permit the weights list to be formed with zero-length weights vectors}
  \item{from_mat2listw}{default FALSE, set TRUE if called from \code{mat2listw}}
}

\value{
\code{listw2sn()}returns a data frame with three columns, and with class \code{spatial.neighbour}:
  \item{from}{region number id for the start of the link (S-PLUS row.id)}
  \item{to}{region number id for the end of the link (S-PLUS col.id)}
  \item{weights}{weight for this link}
%\code{logSpwdet} returns log det(I - rho * W).
}


\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{nb2listw}}}

\examples{
columbus <- st_read(system.file("shapes/columbus.shp", package="spData")[1], quiet=TRUE)
col.gal.nb <- read.gal(system.file("weights/columbus.gal", package="spData")[1])
col.listw <- nb2listw(col.gal.nb)
col.listw$neighbours[[1]]
col.listw$weights[[1]]
col.sn <- listw2sn(col.listw)
str(col.sn)
}
\keyword{spatial}
