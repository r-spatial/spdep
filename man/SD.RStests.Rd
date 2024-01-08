\name{SD.RStests}
\alias{SD.RStests}

\title{
Rao's score and adjusted Rao's score tests of linear hypotheses for spatial Durbin and spatial Durbin error models
}
\description{
Rao's score and adjusted Rao's score tests of linear hypotheses applied to a fitted linear model to examine whether either the spatially lagged dependent variable \code{lag} or the spatially lagged independent variable(s) \code{WX} should be included in the model, or both (SDM). Adjusted tests are provided for \code{lag} and \code{WX} adapting to the presence of the other, and a joint test for both. The joint test is equal to the unadjusted of one plus the adjusted of the other. In addition, draft tests are added from Koley (2024, section 6) for spatial Durbin error models to examine whether either the spatially lagged error \code{err} or the spatially lagged independent variable(s) \code{WX} should be included in the model, or both (SDEM); because of orthogonality, no adjusted tests are required.
}
\usage{
SD.RStests(model, listw, zero.policy = attr(listw, "zero.policy"), test = "SDM",
 Durbin = TRUE)
}

\arguments{
  \item{model}{an object of class \code{lm} returned by \code{lm}}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw},
expected to be row-standardised (W-style)}
  \item{zero.policy}{default \code{attr(listw, "zero.policy")} as set when \code{listw} was created, if attribute not set, use global option value; if TRUE assign zero to the lagged value of zones without
neighbours, if FALSE assign NA}
  \item{test}{test=\dQuote{SDM} computes the SDM tests, a character vector of tests requested chosen from SDM_RSlag, SDM_adjRSlag, SDM_RSWX, SDM_adjRSWX, SDM_Joint, test=\dQuote{SDEM} computes the SDEM tests, a character vector of tests requested chosen from SDEM_RSerr, SDEM_RSWX, SDEM_Joint; test=\dQuote{all} computes all the tests}
  \item{Durbin}{default TRUE for Durbin models including WX; if TRUE, full spatial Durbin model; if a formula object, the subset of explanatory variables to lag}
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
Malabika Koley and Anil K. Bera (2024) To use, or not to use the spatial Durbin model? – that is the question, Spatial Economic Analysis, 19:1, 30-56, \doi{10.1080/17421772.2023.2256810}; Malabika Koley (2024) Specification Testing under General Nesting Spatial Model, \url{https://sites.google.com/view/malabikakoley/research}.
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}, Malabika Koley and Anil K. Bera}
\note{
The results in the example below agree with those in Table 3, p. 22 in Koley and Bera (2024).
}

\seealso{
\code{\link{lm}}, \code{\link{lm.RStests}}
}
\examples{
columbus <- sf::st_read(system.file("shapes/columbus.shp", package="spData")[1])
col.gal.nb <- read.gal(system.file("weights/columbus.gal", package="spData")[1])
col.listw <- nb2listw(col.gal.nb, style="W")
lm_obj <- lm(CRIME ~ INC + HOVAL, data=columbus)
summary(lm.RStests(lm_obj, col.listw, test="all"))
res <- SD.RStests(lm_obj, col.listw, test="SDM")
summary(res)
all.equal(unname(res$SDM_Joint$statistic),
 unname(res$SDM_RSlag$statistic + res$SDM_adjRSWX$statistic))
all.equal(unname(res$SDM_Joint$statistic),
 unname(res$SDM_adjRSlag$statistic + res$SDM_RSWX$statistic))
res <- SD.RStests(lm_obj, col.listw, test="SDEM")
summary(res)
all.equal(unname(res$SDEM_Joint$statistic),
 unname(res$SDEM_RSerr$statistic + res$SDEM_RSWX$statistic))
summary(SD.RStests(lm_obj, nb2listw(col.gal.nb, style="C"), test="all"))
summary(SD.RStests(lm_obj, col.listw, test="all", Durbin= ~ INC))
lm_obj0 <- lm(I(scale(CRIME)) ~ 0 + I(scale(INC)) + I(scale(HOVAL)),
 data=columbus)
summary(SD.RStests(lm_obj0, col.listw, test="all"))
columbusNA <- columbus
columbusNA$HOVAL[15] <- NA
lm_objNA <- lm(CRIME ~ INC + HOVAL, data=columbusNA)
summary(SD.RStests(lm_objNA, col.listw, test="all"))
}
\keyword{spatial}

