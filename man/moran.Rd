% Copyright 2001 by Roger S. Bivand
\name{moran}
\alias{moran}
\title{Compute Moran's I}
\description{
  A simple function to compute Moran's I, called by \code{moran.test} and \code{moran.mc};
\deqn{I = \frac{n}{\sum_{i=1}^{n}\sum_{j=1}^{n}w_{ij}}
\frac{\sum_{i=1}^{n}\sum_{j=1}^{n}w_{ij}(x_i-\bar{x})(x_j-\bar{x})}{\sum_{i=1}^{n}(x_i - \bar{x})^2}
}{I = (n sum_i sum_j w_ij (x_i - xbar) (x_j - xbar)) / (S0 sum_i (x_i - xbar)^2)}
}
\usage{
moran(x, listw, n, S0, zero.policy=NULL, NAOK=FALSE)
}
\arguments{
  \item{x}{a numeric vector the same length as the neighbours list in listw}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{n}{number of zones}
  \item{S0}{global sum of weights}
  \item{zero.policy}{default NULL, use global option value; if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA}
  \item{NAOK}{if 'TRUE' then any 'NA' or 'NaN' or 'Inf' values in x are passed on to the foreign function.  If 'FALSE', the presence of 'NA' or 'NaN' or 'Inf' values is regarded as an error.}
}

\value{
a list of
  \item{I}{Moran's I}
  \item{K}{sample kurtosis of x}
}
\references{Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 17.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{moran.test}}, \code{\link{moran.mc}}}

\examples{
data(oldcol)
col.W <- nb2listw(COL.nb, style="W")
crime <- COL.OLD$CRIME
str(moran(crime, col.W, length(COL.nb), Szero(col.W)))
is.na(crime) <- sample(1:length(crime), 10)
str(moran(crime, col.W, length(COL.nb), Szero(col.W), NAOK=TRUE))
}
\keyword{spatial}
