% Copyright 2004 by Roger S. Bivand
\name{nb2blocknb}
\alias{nb2blocknb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Block up neighbour list for location-less observations}
\description{
  The function blocks up a neighbour list for known spatial locations to
create a new neighbour list for multiple location-less observations know to
belong to the spatial locations, using the identification tags of the
locations as the key. 
}
\usage{
nb2blocknb(nb=NULL, ID, row.names = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nb}{an object of class \code{nb} with a list of integer vectors containing neighbour region number ids; if null, an nb object with no neighbours is created the length of \code{unique(as.character(ID))}}
  \item{ID}{identification tags of the locations for the location-less observations;

 \code{sort(unique(as.character(ID)))} must be identical to 

\code{sort(as.character(attr(nb, "region.id")))};

 same length as \code{row.names} if provided.}
  \item{row.names}{character vector of observation ids to be added to the neighbours list as attribute \code{region.id}, default \code{seq(1, nrow(x))}; same length as \code{ID} if provided.}
}
\details{
Assume that there is a list of unique
locations, then a neighbour list can build for that, to create an input
neighbour list. This needs to be "unfolded", so that observations
belonging to each unique location are observation neighbours, and
observations belonging to the location neighbours of the unique location
in question are also observation neighbours, finally removing the
observation itself (because it should not be its own neighbour). This
scenario also arises when say only post codes are available, and some post
codes contain multiple observations, where all that is known is that they belong
to a specific post code, not where they are located within it (given that 
the post code locations are known).
}

\value{
The function returns an object of class \code{nb} with a list of integer vectors containing neighbour observation number ids.
}

\author{Roger Bivand \email{Roger.Bivand@nhh.no}}


\seealso{\code{\link{knn2nb}}, \code{\link{dnearneigh}}, \code{\link{cell2nb}}, \code{\link{tri2nb}}, \code{\link{poly2nb}}}

\examples{
\dontrun{
data(boston, package="spData")
summary(as.vector(table(boston.c$TOWN)))
townaggr <- aggregate(boston.utm, list(town=boston.c$TOWN), mean)
block.rel <- graph2nb(relativeneigh(as.matrix(townaggr[,2:3])),
 as.character(townaggr[,1]), sym=TRUE)
block.rel
print(is.symmetric.nb(block.rel))
plot(block.rel, as.matrix(townaggr[,2:3]))
points(boston.utm, pch=18, col="lightgreen")
block.nb <- nb2blocknb(block.rel, as.character(boston.c$TOWN))
block.nb
print(is.symmetric.nb(block.nb))
plot(block.nb, boston.utm)
points(boston.utm, pch=18, col="lightgreen")
n.comp.nb(block.nb)$nc
moran.test(boston.c$CMEDV, nb2listw(boston.soi))
moran.test(boston.c$CMEDV, nb2listw(block.nb))
block.nb <- nb2blocknb(NULL, as.character(boston.c$TOWN))
block.nb
print(is.symmetric.nb(block.nb))
plot(block.nb, boston.utm)
n.comp.nb(block.nb)$nc
moran.test(boston.c$CMEDV, nb2listw(block.nb, zero.policy=TRUE), zero.policy=TRUE)
}
}
\keyword{spatial}

