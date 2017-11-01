\name{nb2INLA}
\alias{nb2INLA}
\title{Output spatial neighbours for INLA}
\description{
  Output spatial neighbours for INLA
}
\usage{
nb2INLA(file, nb)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{file where adjacency matrix will be stored}
  \item{nb}{an object of class \code{nb}}
}

\value{
  Nothing is returned but a file will be created with the representation of
the adjacency matrix as required by INLA for its spatial models.
}
\references{http://www.r-inla.org}
\author{Virgilio Gomez-Rubio}

%\seealso{\code{\link[INLA]{inla}}}
\examples{
if (require(rgdal, quietly=TRUE)) {
example(columbus, package="spData")
td <- tempdir()
x <- nb2INLA(paste(td, "columbus-INLA.adj", sep="/"), col.gal.nb)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}

