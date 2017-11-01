% Copyright 2004 by Roger S. Bivand
\name{subset.listw}
\alias{subset.listw}
\title{Subset a spatial weights list}
\description{
  The function subsets a spatial weights list, retaining objects for which the subset argument vector is TRUE. At present it will only subset non-general weights lists (that is those created by \code{nb2listw} with \code{glist=NULL}).}
\usage{
\method{subset}{listw}(x, subset, zero.policy = NULL, ...)
}

\arguments{
  \item{x}{an object of class \code{listw}}
  \item{subset}{logical expression}
  \item{zero.policy}{default NULL, use global option value; if FALSE stop with error for any empty neighbour sets, if TRUE permit the weights list to be formed with zero-length weights vectors - passed through to \code{nb2listw}}
  \item{\dots}{generic function pass-through}
}
\value{
  The function returns an object of class \code{listw} with component \code{style} the same as the input object, component \code{neighbours} a list of integer vectors containing neighbour region number ids (compacted to run from 1:number of regions in subset), and component \code{weights} as the weights computed for \code{neighbours} using \code{style}.
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}
\seealso{\code{\link{nb2listw}}, \code{\link{subset.nb}}}


\examples{
if (require(rgdal, quietly=TRUE)) {
example(columbus, package="spData")
to.be.dropped <- c(31, 34, 36, 39, 42, 46)
pre <- nb2listw(col.gal.nb)
print(pre)
post <- subset(pre, !(1:length(col.gal.nb) \%in\% to.be.dropped))
print(post)
}}
\keyword{spatial}
