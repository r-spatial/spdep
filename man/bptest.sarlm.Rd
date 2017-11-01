% Copyright 2004-2012 by Roger S. Bivand
\encoding{UTF-8}
\name{bptest.sarlm}
\alias{bptest.sarlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Breusch-Pagan test for spatial models}
\description{
  Performs the Breusch-Pagan test for heteroskedasticity on the least squares fit of the spatial models taking the spatial coefficients rho or lambda into account. This function is a copy of the \code{bptest} function in package "lmtest", modified to use objects returned by spatial simultaneous autoregressive models.
}
\usage{
bptest.sarlm(object, varformula=NULL, studentize = TRUE, data=list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{"sarlm"} from \code{errorsarlm()} or \code{lagsarlm()}.}
  \item{varformula}{a formula describing only the potential explanatory variables for the variance (no dependent variable needed). By default the same explanatory variables are taken as in the main regression model}
  \item{studentize}{logical. If set to \code{TRUE} Koenker's studentized
   version of the test statistic will be used.}
  \item{data}{an optional data frame containing the variables in the varformula}
}
\details{
   Asymptotically this corresponds to the test given by Anselin (1988), but is not exactly the same. The studentized version is more conservative and perhaps to be prefered. The residuals, and for spatial error models the RHS variables, are adjusted for the spatial coefficient, as suggested bt Luc Anselin (personal communication).

It is also technically possible to make heteroskedasticity corrections to standard error estimates by using the \dQuote{lm.target} component of \code{sarlm} objects - using functions in the lmtest and sandwich packages.
}

\value{
 A list with class \code{"htest"} containing the following components:
  \item{statistic}{the value of the test statistic.}
  \item{p.value}{the p-value of the test.}
  \item{parameter}{degrees of freedom (wrongly reported if varformula given before 0.5-44.}
  \item{method}{a character string indicating what type of test was
    performed.}
}
\references{T.S. Breusch & A.R. Pagan (1979),
A Simple Test for Heteroscedasticity and Random Coefficient Variation.
\emph{Econometrica} \bold{47}, 1287--1294

W. Krämer & H. Sonnberger (1986),
\emph{The Linear Regression Model under Test}. Heidelberg: Physica.

L. Anselin (1988) \emph{Spatial econometrics: methods and models.}
Dordrecht: Kluwer, pp. 121--122.
}
\author{Torsten Hothorn \email{Torsten.Hothorn@rzmail.uni-erlangen.de} and
 Achim Zeileis \email{zeileis@ci.tuwien.ac.at}, modified by Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{errorsarlm}}, \code{\link{lagsarlm}}}
\examples{
example(columbus)
error.col <- errorsarlm(CRIME ~ HOVAL + INC, data=columbus,
 nb2listw(col.gal.nb))
bptest.sarlm(error.col)
bptest.sarlm(error.col, studentize=FALSE)
\dontrun{
lm.target <- lm(error.col$tary ~ error.col$tarX - 1)
if (require(lmtest) && require(sandwich)) {
  coeftest(lm.target, vcov=vcovHC(lm.target, type="HC0"), df=Inf)
}
}
}
\keyword{spatial}

