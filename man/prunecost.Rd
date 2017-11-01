\name{prunecost}
\alias{prunecost}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute cost of prune each edge}
\description{
  If any edge are dropped, the MST are pruned. This generate a two
  subgraphs. So, it makes a tree graphs and tree dissimilarity values
  are computed, one for each graph. The dissimilarity is the sum over
  sqared differences between the observactions in the nodes and mean
  vector of observations in the graph. The dissimilarity of original
  graph and the sum of dissimilarity of subgraphs are returned.
}
\usage{
prunecost(edges, data, method = c("euclidean", "maximum", "manhattan", 
    "canberra", "binary", "minkowski", "mahalanobis"), 
    p = 2, cov, inverted = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{edges}{A matrix with 2 colums with each row is one edge}
  \item{data}{A data.frame with observations in the nodes.}
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
%%\details{}
\value{
  A vector with the differences between the dissimilarity of all nodes
  and the dissimilarity sum of all subgraphs obtained by pruning one
  edge each time.
}
%%\references{ ~put references to the literature/web site here ~ }
\author{Elias T. Krainski and Renato M. Assuncao}
%%\note{}
\seealso{See Also as \code{\link{prunemst}}}
\examples{
d <- data.frame(a=-2:2, b=runif(5))
e <- matrix(c(1,2, 2,3, 3,4, 4,5), ncol=2, byrow=TRUE)

sum(sweep(d, 2, colMeans(d))^2)

prunecost(e, d)
}
\keyword{graphs}
\keyword{cluster}% __ONLY ONE__ keyword per line
