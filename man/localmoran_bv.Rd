\name{localmoran_bv}
\alias{localmoran_bv}
\title{Compute the Local Bivariate Moran's I Statistic}
\usage{
localmoran_bv(x, y, listw, nsim = 199, scale = TRUE, alternative="two.sided",
 iseed=1L, no_repeat_in_row=FALSE, zero.policy=attr(listw, "zero.policy"))
}
\arguments{
\item{x}{a numeric vector of same length as \code{y}.}

\item{y}{a numeric vector of same length as \code{x}.}

\item{listw}{a listw object for example as created by \code{nb2listw()}.}

\item{nsim}{the number of simulations to run.}

\item{scale}{default \code{TRUE}.}

\item{alternative}{a character string specifying the alternative hypothesis, must be one of "greater" (default), "two.sided", or "less".}
\item{iseed}{default NULL, used to set the seed for possible parallel RNGs.}
\item{no_repeat_in_row}{default \code{FALSE}, if \code{TRUE}, sample conditionally in each row without replacements to avoid duplicate values, \url{https://github.com/r-spatial/spdep/issues/124}}
  \item{zero.policy}{default default \code{attr(listw, "zero.policy")} as set when \code{listw} was created, if attribute not set, use global option value; if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA}
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
# load columbus data
columbus <- st_read(system.file("shapes/columbus.shp", package="spData"))
nb <- poly2nb(columbus)
listw <- nb2listw(nb)
set.seed(1)
(res <- localmoran_bv(columbus$CRIME, columbus$INC, listw, nsim = 499))
columbus$hs <- hotspot(res, Prname="Pr(folded) Sim", cutoff=0.05,
 quadrant.type="pysal", p.adjust="none")
if (require("tmap", quietly=TRUE)) {
  tm_shape(columbus) + tm_fill("hs")
}
moran.plot(x=columbus$CRIME, y=columbus$INC, listw=listw)
}
\references{
Anselin, Luc, Ibnu Syabri, and Oleg Smirnov. 2002. “Visualizing Multivariate Spatial Correlation with Dynamically Linked Windows.” In New Tools for Spatial Data Analysis: Proceedings of the Specialist Meeting, edited by Luc Anselin and Sergio Rey. University of California, Santa Barbara: Center for Spatially Integrated Social Science (CSISS).
}

\author{Josiah Parry \email{josiah.parry@gmail.com}}
