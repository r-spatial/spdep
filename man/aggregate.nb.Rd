\name{aggregate.nb}
\alias{aggregate.nb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Aggregate a spatial neighbours object}
\description{
  The method aggregates a spatial neighbours object, creating a new object listing the neighbours of the aggregates.
}
\usage{
\method{aggregate}{nb}(x, IDs, remove.self = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an nb neighbour object}
  \item{IDs}{a character vector of IDs grouping the members of the neighbour object}
  \item{remove.self}{default TRUE: remove self-neighbours resulting from aggregation}
  \item{\dots}{unused - arguments passed through}
}

\value{
  an nb neighbour object
}

\author{Roger Bivand \email{Roger.Bivand@nhh.no}}
\note{Method suggested by Roberto Patuelli}

\examples{
data(used.cars, package="spData")
data(state)
cont_st <- match(attr(usa48.nb, "region.id"), state.abb)
cents <- as.matrix(as.data.frame(state.center))[cont_st,]
opar <- par(mfrow=c(2,1))
plot(usa48.nb, cents, xlim=c(-125, -65), ylim=c(25, 50))
IDs <- as.character(state.division[cont_st])
agg_cents <- aggregate(cents, list(IDs), mean)
agg_nb <- aggregate(usa48.nb, IDs)
plot(agg_nb, agg_cents[, 2:3], xlim=c(-125, -65), ylim=c(25, 50))
text(agg_cents[, 2:3], agg_cents[, 1], cex=0.6)
par(opar)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}

