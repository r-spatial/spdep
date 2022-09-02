\name{moran_bv}
\alias{moran_bv}
\title{Compute the Global Bivariate Moran's I}
\usage{
moran_bv(x, y, listw, nsim = 99, scale = TRUE)
}
\arguments{
\item{x}{a numeric vector of same length as \code{y}.}

\item{y}{a numeric vector of same length as \code{x}.}

\item{listw}{a listw object for example as created by \code{nb2listw()}.}

\item{nsim}{the number of simulations to run.}

\item{scale}{default \code{TRUE}.}
}
\value{
a named list with two elements \code{Ib} and \code{p_sim} containing the bivariate Moran's I and simulated p-value respectively.
}
\description{
Given two continuous numeric variables, calculate the bivariate Moran's I. See details for more.
}
\details{
The Global Bivariate Moran is defined as

\eqn{
I_B = \frac{\Sigma_i(\Sigma_j{w_{ij}y_j\times x_i})}{\Sigma_i{x_i^2}}
}

It is important to note that this is a measure of autocorrelation of X
with the spatial lag of Y. As such, the resultant measure may overestimate the amount of
spatial autocorrelation which may be a product of the inherent correlation of X and Y.
}
\examples{
data(boston, package = "spData")
x <- boston.c$CRIM
y <- boston.c$NOX
listw <- nb2listw(boston.soi)
moran_bv(x, y, listw)
}
\references{
\href{https://geodacenter.github.io/workbook/5b_global_adv/lab5b.html}{Global Spatial Autocorrelation (2): Bivariate, Differential and EB Rate Moran Scatter Plot, Luc Anselin}
}

\author{Josiah Parry \email{josiah.parry@gmail.com}}
