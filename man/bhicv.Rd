\name{bhicv}

\alias{bhicv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Data set with 4 life condition indices of Belo Horizonte region}
\description{
  The data are collected inthe  Atlas of condition indices published 
  by the Joao Pinheiro Foundation and UNDP.
}

\format{
  A shape polygon object with seven variables:
  \describe{
  \item{id}{The identificator}
  \item{Name}{Name of city}
  \item{Population}{The population of city}
  \item{HLCI}{Health Life Condition Index}
  \item{ELCI}{Education Life Condition Index}
  \item{CLCI}{Children Life Condition Index}
  \item{ELCI}{Economic Life Condition Index}
  }
}
%- maybe also 'usage' for other objects documented here.
%%\format{}
%%\details{}
\examples{
if (require(rgdal, quietly=TRUE)) {
bh <- readOGR(system.file("etc/shapes/bhicv.shp",
      package="spdep")[1])
}
}
\keyword{data}% at least one, from doc/KEYWORDS
