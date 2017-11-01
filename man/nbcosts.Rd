\name{nbcosts}
\alias{nbcosts}
\alias{nbcost}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute cost of edges}
\description{
  The cost of each edge is the distance between it nodes. This function
  compute this distance using a data.frame with observations vector in
  each node.
}
\usage{
nbcost(data, id, id.neigh,  method = c("euclidean", "maximum", 
    "manhattan", "canberra", "binary", "minkowski", "mahalanobis"),
    p = 2, cov, inverted = FALSE)
nbcosts(nb, data,  method = c("euclidean", "maximum", 
    "manhattan", "canberra", "binary", "minkowski", "mahalanobis"),
    p = 2, cov, inverted = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nb}{An object of \code{nb} class. See \code{\link{poly2nb}} for
    details.} 
  \item{data}{A matrix with observations in the nodes.}
  \item{id}{Node index to compute the cost}
  \item{id.neigh}{Idex of neighbours nodes of node \code{id}}
  \item{method}{Character or function to declare distance method.
    If \code{method} is character, method must be "mahalanobis" or
    "euclidean", "maximum", "manhattan", "canberra", "binary"
    or "minkowisk".
    If \code{method} is one of "euclidean", "maximum",
    "manhattan", "canberra", "binary" or "minkowisk", see
    \code{\link{dist}} for details,
     because this function as used to compute the distance.
     If \code{method="mahalanobis"}, the mahalanobis distance
     is computed between neighbour areas.
     If \code{method} is a \code{function}, this function is
     used to compute the distance.}
  \item{p}{The power of the Minkowski distance.}
  \item{cov}{The covariance matrix used to compute the mahalanobis 
   distance.}
  \item{inverted}{logical.  If 'TRUE', 'cov' is supposed to contain the
          inverse of the covariance matrix.}
}
%\details{}
\value{
  A object of \code{nbdist} class. See \code{\link{nbdists}} for
  details.
}
\note{The neighbours must be a connected graph.}
%%\references{ ~put references to the literature/web site here ~ }
\author{Elias T. Krainski and Renato M. Assuncao}
%%\note{}
\seealso{See Also as \code{\link{nbdists}}, \code{\link{nb2listw}}}
%%\examples{}
\keyword{cluster}
\keyword{spatial}% __ONLY ONE__ keyword per line
