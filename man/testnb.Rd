% Copyright 2001 by Roger S. Bivand
\name{is.symmetric.nb}
\alias{is.symmetric.nb}
\alias{sym.attr.nb}
\alias{make.sym.nb}
\alias{old.make.sym.nb}
\alias{is.symmetric.glist}

\title{Test a neighbours list for symmetry}
\description{
Checks a neighbours list for symmetry/transitivity (if i is a neighbour of j,
then j is a neighbour of i). This holds for distance and contiguity based
neighbours, but not for k-nearest neighbours. The helper function
\code{sym.attr.nb()} calls \code{is.symmetric.nb()} to set the \code{sym}
attribute if needed, and \code{make.sym.nb} makes a non-symmetric list symmetric by adding neighbors. \code{is.symmetric.glist} checks a list of general weights corresponding to neighbours for symmetry for symmetric neighbours.
}
\usage{
is.symmetric.nb(nb, verbose = NULL, force = FALSE)
sym.attr.nb(nb)
make.sym.nb(nb)
old.make.sym.nb(nb)
is.symmetric.glist(nb, glist)
}

\arguments{
  \item{nb}{an object of class \code{nb} with a list of integer vectors containing neighbour region number ids.}
  \item{verbose}{default NULL, use global option value; if TRUE prints non-matching pairs}
  \item{force}{do not respect a neighbours list \code{sym} attribute and test anyway}
  \item{glist}{list of general weights corresponding to neighbours}
}
\value{
TRUE if symmetric, FALSE if not; is.symmetric.glist returns a value with an attribute, "d", indicating for failed symmetry the largest failing value.
}

\note{
A new version of \code{make.sym.nb} by Bjarke Christensen is now included. The older version has been renamed \code{old.make.sym.nb}, and their comparison constitutes a nice demonstration of vectorising speedup using \code{sapply} and \code{lapply} rather than loops. When any no-neighbour observations are present, \code{old.make.sym.nb} is used.
}

\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{read.gal}}}

\examples{
if (require(rgdal, quietly=TRUE)) {
example(columbus, package="spData")
coords <- coordinates(columbus)
ind <- sapply(slot(columbus, "polygons"), function(x) slot(x, "ID"))
print(is.symmetric.nb(col.gal.nb, verbose=TRUE, force=TRUE))
k4 <- knn2nb(knearneigh(coords, k=4), row.names=ind)
k4 <- sym.attr.nb(k4)
print(is.symmetric.nb(k4))
k4.sym <- make.sym.nb(k4)
print(is.symmetric.nb(k4.sym))
}}
\keyword{spatial}
