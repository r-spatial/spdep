\name{nb2WB}
\alias{nb2WB}
\alias{listw2WB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Output spatial weights for WinBUGS}
\description{
  Output spatial weights for WinBUGS
}
\usage{
nb2WB(nb)
listw2WB(listw) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nb}{an object of class \code{nb}}
  \item{listw}{a \code{listw} object from for example \code{nb2listw}}
}

\value{
  A list suitable for convering using \code{dput} for WinBUGS
}
\references{http://www.mrc-bsu.cam.ac.uk/bugs/winbugs/geobugs12manual.pdf}
\author{Virgilio Gomez-Rubio}

\seealso{\code{\link[base]{dput}}}
\examples{
if (require(rgdal, quietly=TRUE)) {
example(columbus, package="spData")
x <- nb2WB(col.gal.nb)
dput(x, control=NULL)
x <- listw2WB(nb2listw(col.gal.nb))
dput(x, control=NULL)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}

