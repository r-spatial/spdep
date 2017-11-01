\name{plot.skater}
\alias{plot.skater}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot the object of skater class}
\description{
  This function displays the results of the skater function. The subgraphs
  are plotted with different colours. 
}
\usage{
\method{plot}{skater}(x, coords, label.areas = NULL, 
   groups.colors, cex.circles = 1, cex.labels = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of \code{skater} class.}
  \item{coords}{A matrix of two colums with coordinates of nodes.}
  \item{label.areas}{A vector of labels of nodes.}
  \item{groups.colors}{A vector with colors of groups ou sub-graphs.}
  \item{cex.circles}{The length of circles with represent the nodes.}
  \item{cex.labels}{The length of labels of nodes.}
  \item{\dots}{Further arguments passed to plotting funcitons.}
}
\author{Elias T. Krainski and Renato M. Assuncao}
\seealso{See Also as \code{\link{skater}} and \code{\link{mstree}}}
\examples{
### see example in the skater function documentation
}
\keyword{hplot}
\keyword{cluster}
