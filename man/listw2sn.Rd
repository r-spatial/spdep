% Copyright 2001-3 by Roger S. Bivand
\name{listw2sn}
\alias{listw2sn}
\alias{sn2listw}
\alias{as.spam.listw}
\title{Spatial neighbour sparse representation}
\description{
The function makes a \code{"spatial neighbour"} object representation
(similar to the S-PLUS spatial statististics module representation of a
\code{"listw"} spatial weights object. \code{sn2listw()} is the inverse function to \code{listw2sn()}, creating a \code{"listw"} object from a \code{"spatial neighbour"} object. The \code{as.spam.listw} method converts a \code{"listw"} object to a sparse matrix as defined in the \pkg{spam} package, using \code{listw2sn()}.
}
\usage{
listw2sn(listw)
sn2listw(sn)
as.spam.listw(listw)
}
\arguments{
  \item{listw}{a \code{listw} object from for example \code{nb2listw}}
  \item{sn}{a \code{spatial.neighbour} object}
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
example(columbus)
col.listw <- nb2listw(col.gal.nb)
col.listw$neighbours[[1]]
col.listw$weights[[1]]
col.sn <- listw2sn(col.listw)
str(col.sn)
\dontrun{
col.sp <- as.spam.listw(col.listw)
str(col.sp)
}
}
\keyword{spatial}
