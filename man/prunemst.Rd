\name{prunemst}
\alias{prunemst}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Prune a Minimun Spanning Tree}
\description{
  This function deletes a first edge and makes two subsets of edges. Each
  subset is a Minimun Spanning Treee.
}
\usage{
prunemst(edges, only.nodes = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{edges}{A matrix with two colums with each row is one edge}
  \item{only.nodes}{If \code{only.nodes=FALSE}, return a edges and nodes
  of each MST resulted. If \code{only.nodes=TRUE}, return a two sets of
  nodes. Defalt is TRUE} 
}
%%\details{}
\value{
  A list of length two. If \code{only.nodes=TRUE} each element is a
  vector of nodes. If \code{only.nodes=FALSE} each element is a list with
  nodes and edges.
}
%%\references{ ~put references to the literature/web site here ~ }
\author{Elias T. Krainski and Renato M. Assuncao}
%%\note{ ~~further notes~~ }
\seealso{ See Also as \code{\link{mstree}} }
\examples{
e <- matrix(c(2,3, 1,2, 3,4, 4,5), ncol=2, byrow=TRUE)
e
prunemst(e)
prunemst(e, only.nodes=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{tree}
\keyword{cluster}% __ONLY ONE__ keyword per line
