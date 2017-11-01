% Copyright 2001 by Roger S. Bivand
\name{airdist}
\alias{airdist}
\title{Measure distance from plot}
\description{
Measure a distance between two points on a plot using \code{locator}; the function checks \code{par("plt")} and \code{par("usr")} to try to ensure that the aspect ratio y/x is 1, that is that the units of measurement in both x and y are equivalent.
}
\usage{
airdist(ann=FALSE)
}
\arguments{
  \item{ann}{annotate the plot with line measured and distance}
}
\value{
a list with members:
  \item{dist}{distance measured}
  \item{coords}{coordinates between which distance is measured}
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{locator}}}

\keyword{spatial}
