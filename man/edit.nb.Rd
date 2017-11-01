% Copyright 2001-2010 by Roger S. Bivand
\name{edit.nb}
\alias{edit.nb}
\title{Interactive editing of neighbours lists}
\description{
The function provides simple interactive editing of neighbours lists to
allow unneeded links to be deleted, and missing links to be inserted. It uses \code{identify} to pick the
endpoints of the link to be deleted or added, and asks for confirmation before
committing. If the result is not assigned to a new object, the editing
will be lost - as in \code{edit}.
}
\usage{
\method{edit}{nb}(name, coords, polys=NULL, ..., use_region.id=FALSE)
}

\arguments{
  \item{name}{an object of class \code{nb}}
  \item{coords}{matrix of region point coordinates; if missing and polys= inherits from \code{SpatialPolygons}, the label points of that object are used}
  \item{polys}{if polygon boundaries supplied, will be used as background; must inherit from \code{SpatialPolygons}}
  \item{...}{further arguments passed to or from other methods}
  \item{use_region.id}{default \code{FALSE}, in \code{identify} use 1-based observation numbers, otherwise use the \code{nb} \code{region.id} attribute values}
}
\value{
The function returns an object of class \code{nb} with the edited list of integer vectors containing neighbour region number ids, with added attributes tallying the added and deleted links.
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{summary.nb}}, \code{\link{plot.nb}}}

\examples{
\dontrun{
if (require(rgdal, quietly=TRUE)) {
example(columbus, package="spData")
class(columbus)
nnb1 <- edit.nb(col.gal.nb, polys=columbus)
}}}

\keyword{spatial}
