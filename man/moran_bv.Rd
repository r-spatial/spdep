\name{moran_bv}
\alias{moran_bv}
\title{Compute the Global Bivariate Moran's I}
\usage{
moran_bv(x, y, listw, nsim = 499, scale = TRUE)
}
\arguments{
\item{x}{a numeric vector of same length as \code{y}.}

\item{y}{a numeric vector of same length as \code{x}.}

\item{listw}{a listw object for example as created by \code{nb2listw()}.}

\item{nsim}{the number of simulations to run.}

\item{scale}{default \code{TRUE}.}
}
\value{
An object of class \code{"boot"}, with the observed statistic in component \code{t0}.
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
spatial autocorrelation which may be a product of the inherent correlation of X and Y. The output object is of class \code{"boot"}, so that plots and confidence intervals are available using appropriate methods.
}
\examples{
data(boston, package = "spData")
x <- boston.c$CRIM
y <- boston.c$NOX
listw <- nb2listw(boston.soi)
set.seed(1)
res_xy <- moran_bv(x, y, listw, nsim=499)
res_xy$t0
boot.ci(res_xy, conf=c(0.99, 0.95, 0.9), type="basic")
plot(res_xy)
set.seed(1)
lee_xy <- lee.mc(x, y, listw, nsim=499, return_boot=TRUE)
lee_xy$t0
boot.ci(lee_xy, conf=c(0.99, 0.95, 0.9), type="basic")
plot(lee_xy)
set.seed(1)
res_yx <- moran_bv(y, x, listw, nsim=499)
res_yx$t0
boot.ci(res_yx, conf=c(0.99, 0.95, 0.9), type="basic")
plot(res_yx)
set.seed(1)
lee_yx <- lee.mc(y, x, listw, nsim=499, return_boot=TRUE)
lee_yx$t0
boot.ci(lee_yx, conf=c(0.99, 0.95, 0.9), type="basic")
plot(lee_yx)
}
\references{
Wartenberg, D. (1985), Multivariate Spatial Correlation: A Method for Exploratory Geographical Analysis. Geographical Analysis, 17: 263-283. \doi{10.1111/j.1538-4632.1985.tb00849.x}
}

\author{Josiah Parry \email{josiah.parry@gmail.com}}
