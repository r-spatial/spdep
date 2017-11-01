\name{plot.mst}
\alias{plot.mst}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot the Minimum Spanning Tree}
\description{
  This function plots a MST, the nodes are circles and the edges are
  segments. 
}
\usage{
\method{plot}{mst}(x, coords, label.areas = NULL, 
   cex.circles = 1, cex.labels = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Object of \code{mst} class.}
  \item{coords}{A two column matrix with the coordinates of nodes.}
  \item{label.areas}{A vector with the labels of nodes}
  \item{cex.circles}{The length of circles to plot.}
  \item{cex.labels}{The length of nodes labels ploted.}
  \item{\dots}{Further arguments passed to plotting funcitons.}
}
\author{Elias T. Krainski and Renato M. Assuncao}
\seealso{See Also as \code{\link{skater}} and \code{\link{mstree}}}
\examples{
### see example in mstree function documentation
}
\keyword{hplot}
\keyword{tree}
