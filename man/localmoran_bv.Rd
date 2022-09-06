\name{localmoran_bv}
\alias{localmoran_bv}
\title{Compute the Local Bivariate Moran's I Statistic}
\usage{
localmoran_bv(x, y, listw, nsim = 199, scale = TRUE)
}
\arguments{
\item{x}{a numeric vector of same length as \code{y}.}

\item{y}{a numeric vector of same length as \code{x}.}

\item{listw}{a listw object for example as created by \code{nb2listw()}.}

\item{nsim}{the number of simulations to run.}

\item{scale}{default \code{TRUE}.}
}
\value{
a \code{data.frame} containing two columns \code{Ib} and \code{p_sim} containing the local bivariate Moran's I and simulated p-values respectively.
}
\description{
Given two continuous numeric variables, calculate the bivariate Local Moran's I.
}
\details{
The Bivariate Local Moran, like its global counterpart, evaluates the value
of x at observation i with its spatial neighbors' value of y. The value of \eqn{I_i^B} is xi * Wyi. Or, in simpler words, the local bivariate Moran is the result of multiplying x by the spatial lag of y. Formally it is defined as

\eqn{
I_i^B= cx_i\Sigma_j{w_{ij}y_j}
}
}
\examples{
# load boston data
data("boston", package = "spData")
x <- boston.c$CRIM
y <- boston.c$NOX
# convert the nb listw to listw
listw <- nb2listw(boston.soi)

localmoran_bv(x, y, listw, nsim = 99)
}
\references{
Anselin, Luc, Ibnu Syabri, and Oleg Smirnov. 2002. “Visualizing Multivariate Spatial Correlation with Dynamically Linked Windows.” In New Tools for Spatial Data Analysis: Proceedings of the Specialist Meeting, edited by Luc Anselin and Sergio Rey. University of California, Santa Barbara: Center for Spatially Integrated Social Science (CSISS).
}

\author{Josiah Parry \email{josiah.parry@gmail.com}}
