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
GDAL37 <- numeric_version(unname(sf::sf_extSoftVersion()["GDAL"]), strict=FALSE)
(GDAL37 <- ifelse(is.na(GDAL37), FALSE, GDAL37 >= "3.7.0"))
file <- "etc/shapes/bhicv.gpkg.zip"
zipfile <- system.file(file, package="spdep")
if (GDAL37) {
    bh <- st_read(zipfile)
} else {
    td <- tempdir()
    bn <- sub(".zip", "", basename(file), fixed=TRUE)
    target <- unzip(zipfile, files=bn, exdir=td)
    bh <- st_read(target)
}
}
\keyword{data}% at least one, from doc/KEYWORDS
