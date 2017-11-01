% Copyright 2001 by Roger S. Bivand
\name{geary}
\alias{geary}
%\alias{geary.intern}
\title{Compute Geary's C}
\description{
  A simple function to compute Geary's C, called by \code{geary.test} and \code{geary.mc};
\deqn{C = \frac{(n-1)}{2\sum_{i=1}^{n}\sum_{j=1}^{n}w_{ij}}
\frac{\sum_{i=1}^{n}\sum_{j=1}^{n}w_{ij}(x_i-x_j)^2}{\sum_{i=1}^{n}(x_i - \bar{x})^2}
}
\code{geary.intern} is an internal function used to vary the similarity
criterion.
}
\usage{
geary(x, listw, n, n1, S0, zero.policy=NULL)
%geary.intern(x, listw, n, zero.policy, type="geary")
}

\arguments{
  \item{x}{a numeric vector the same length as the neighbours list in listw}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{n}{number of zones}
  \item{n1}{n - 1}
  \item{S0}{global sum of weights}
  \item{zero.policy}{default NULL, use global option value; if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA}
%  \item{type}{"geary" or "sokal" criteria for similarity}
}

\value{
a list with
  \item{C}{Geary's C}
  \item{K}{sample kurtosis of x}
}
\references{Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 17.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{geary.test}}, \code{\link{geary.mc}}, 
\code{\link{sp.mantel.mc}}}

\examples{
data(oldcol)
col.W <- nb2listw(COL.nb, style="W")
str(geary(COL.OLD$CRIME, col.W, length(COL.nb), length(COL.nb)-1,
 Szero(col.W)))
}
\keyword{spatial}
