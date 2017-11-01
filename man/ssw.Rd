\name{ssw}
\alias{ssw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute the sum of dissimilarity}
\description{
  This function computes the sum of dissimilarity between each 
  observation and the mean (scalar of vector) of the observations. 
}
\usage{
ssw(data, id, method = c("euclidean", "maximum", 
    "manhattan", "canberra", "binary", "minkowski",
    "mahalanobis"), p = 2, cov, inverted = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A matrix with observations in the nodes.}
  \item{id}{Node index to compute the cost}
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
  A numeric, the sum of dissimilarity between the observations 
  \code{id} of \code{data} and the mean (scalar of vector) of 
  this observations.
}
%%\references{ ~put references to the literature/web site here ~ }
\author{Elias T. Krainski and Renato M. Assuncao}
%%\note{}
\seealso{See Also as \code{\link{nbcost}}}
\examples{
data(USArrests)
n <- nrow(USArrests)
ssw(USArrests, 1:n)
ssw(USArrests, 1:(n/2))
ssw(USArrests, (n/2+1):n)
ssw(USArrests, 1:(n/2)) + ssw(USArrests, (n/2+1):n)
}
\keyword{cluster}
\keyword{multivariate}% __ONLY ONE__ keyword per line
