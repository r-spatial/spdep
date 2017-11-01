% Copyright 2001 by Roger S. Bivand
\name{card}
\alias{card}
\title{Cardinalities for neighbours lists}
\description{
The function tallies the numbers of neighbours of regions in the neighbours
list.
}
\usage{
card(nb)
}
\arguments{
  \item{nb}{a neighbours list object of class \code{nb}}
}
\value{
An integer vector of the numbers of neighbours of regions in the neighbours
list.
}

\details{\dQuote{nb} objects are stored as lists of integer vectors, where the vectors contain either the indices in the range \code{1:n} for \code{n} as {length(nb)} of the neighbours of region \code{i}, or \code{as.integer(0)} to signal no neighbours. The function \code{card(nb)} is used to extract the numbers of neighbours from the \dQuote{nb} object.}


\references{Bivand R, Pebesma EJ, Gomez-Rubio V, (2008) \emph{Applied Spatial Data Analysis with R}, Springer, New York, pp. 239-251; Bivand R, Portnov B, (2004) Exploring spatial data analysis techniques using R: the case of observations with no neighbours. In: Anselin L, Florax R, Rey S, (eds.), \emph{Advances in Spatial Econometrics, Methodology, Tools and Applications}. Berlin: Springer-Verlag, pp. 121-142.}


\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{summary.nb}}}

\examples{
if (require(rgdal, quietly=TRUE)) {
example(columbus, package="spData")
table(card(col.gal.nb))
}}
\keyword{spatial}

