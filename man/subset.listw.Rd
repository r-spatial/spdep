% Copyright 2004 by Roger S. Bivand
\name{subset.listw}
\alias{subset.listw}
\title{Subset a spatial weights list}
\description{
  The function subsets a spatial weights list, retaining objects for which the subset argument vector is TRUE. At present it will only subset non-general weights lists (that is those created by \code{nb2listw} with \code{glist=NULL}).}
\usage{
\method{subset}{listw}(x, subset, zero.policy = attr(x, "zero.policy"), ...)
}

\arguments{
  \item{x}{an object of class \code{listw}}
  \item{subset}{logical expression}
  \item{zero.policy}{default \code{attr(x, "zero.policy")} as set when \code{x} was created, if attribute not set, use global option value; if FALSE stop with error for any empty neighbour sets, if TRUE permit the weights list to be formed with zero-length weights vectors - passed through to \code{nb2listw}}
  \item{\dots}{generic function pass-through}
}
\value{
  The function returns an object of class \code{listw} with component \code{style} the same as the input object, component \code{neighbours} a list of integer vectors containing neighbour region number ids (compacted to run from 1:number of regions in subset), and component \code{weights} as the weights computed for \code{neighbours} using \code{style}. If no-neighbour observations are created by subsetting and \code{zero.policy} in the input weights object was FALSE, it will be set to TRUE and a warning issued.
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}
\seealso{\code{\link{nb2listw}}, \code{\link{subset.nb}}}


\examples{
col.gal.nb <- read.gal(system.file("weights/columbus.gal", package="spData")[1])
to.be.dropped <- c(31, 34, 36, 39, 42, 46)
pre <- nb2listw(col.gal.nb)
print(pre)
post <- subset(pre, !(1:length(col.gal.nb) \%in\% to.be.dropped))
print(post)
columbus <- st_read(system.file("shapes/columbus.shp", package="spData")[1], quiet=TRUE)
nb <- poly2nb(columbus)
lw <- nb2listw(nb, style="W")
attr(lw, "zero.policy")
(lwa <- subset(lw, 1:nrow(columbus) != c(21)))
attr(lwa, "zero.policy")
(lwb <- subset(lw, !(1:nrow(columbus) \%in\% c(21, 36, 39))))
attr(lwb, "zero.policy")
}
\keyword{spatial}
