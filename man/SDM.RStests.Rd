\name{SDM.RStests}
\alias{SDM.RStests}

\title{

}
\description{

}
\usage{
SDM.RStests(model, listw, zero.policy = attr(listw, "zero.policy"), test = "all", Durbin = TRUE)
}

\arguments{
  \item{model}{an object of class \code{lm} returned by \code{lm}}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw},
expected to be row-standardised (W-style)}
  \item{zero.policy}{default \code{attr(listw, "zero.policy")} as set when \code{listw} was created, if attribute not set, use global option value; if TRUE assign zero to the lagged value of zones without
neighbours, if FALSE assign NA}
  \item{test}{test="all" computes all the tests, a character vector of tests requested chosen from SDM_RSlag, SDM_adjRSlag, SDM_RSWX, SDM_adjRSWX, SDM_Joint}
  \item{Durbin}{default TRUE for Durbin models including WX; if TRUE, full spatial Durbin model; if a formula object, the subset of explanatory variables to lag}
}
\details{

}
\value{
A list of class \code{LMtestlist} of \code{htest} objects, each with:
  \item{statistic}{the value of the Lagrange Multiplier test.}
  \item{parameter}{number of degrees of freedom}
  \item{p.value}{the p-value of the test.}
  \item{method}{a character string giving the method used.}
  \item{data.name}{a character string giving the name(s) of the data.}
}
\references{

}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}, Malabika Koley and Anil K. Bera}
\note{

}

\seealso{
\code{\link{lm}}, \code{\link{lm.LMtests}}
}
\examples{
columbus <- sf::st_read(system.file("shapes/columbus.shp", package="spData")[1])
col.gal.nb <- read.gal(system.file("weights/columbus.gal", package="spData")[1])
col.listw <- nb2listw(col.gal.nb, style="W")
lm_obj <- lm(CRIME ~ INC + HOVAL, data=columbus)
summary(lm.LMtests(lm_obj, col.listw, test="all"))
summary(SDM.RStests(lm_obj, col.listw, test="all"))
}
\keyword{spatial}

