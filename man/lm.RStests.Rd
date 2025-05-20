% Copyright 2001-2 by Roger S. Bivand
\name{lm.LMtests}
\alias{lm.LMtests}
\alias{lm.RStests}
\alias{print.RStestlist}
\alias{summary.RStestlist}
\alias{print.RStestlist.summary}
%\alias{tracew}
\title{Rao's score (a.k.a Lagrange Multiplier) diagnostics for spatial dependence in linear models}
\description{
The function reports the estimates of tests chosen among five statistics for
testing for spatial dependence in linear models. The statistics are
the simple RS test for error dependence (\dQuote{RSerr}), the simple RS test
for a missing spatially lagged dependent variable (\dQuote{RSlag}), variants
of these adjusted for the presence of the other (\dQuote{adjRSerr}
tests for error dependence in the possible presence of a missing lagged
dependent variable, \dQuote{adjRSlag} the other way round), and a portmanteau test
(\dQuote{SARMA}, in fact \dQuote{RSerr} + \dQuote{adjRSlag}). Note: from spdep 1.3-2, the tests are re-named \dQuote{RS} - Rao's score tests, rather than \dQuote{LM} - Lagrange multiplier tests to match the naming of tests from the same family in \code{SDM.RStests}.
}
\usage{
lm.RStests(model, listw, zero.policy=attr(listw, "zero.policy"), test="RSerr",
 spChk=NULL, naSubset=TRUE)
lm.LMtests(model, listw, zero.policy=attr(listw, "zero.policy"), test="LMerr",
 spChk=NULL, naSubset=TRUE)
\method{print}{RStestlist}(x, ...)
\method{summary}{RStestlist}(object, p.adjust.method="none", ...)
\method{print}{RStestlist.summary}(x, digits=max(3, getOption("digits") - 2), ...)
%tracew(listw)
}

\arguments{
  \item{model}{an object of class \code{lm} returned by \code{lm}, or optionally a vector of externally calculated residuals (run though \code{na.omit} if any NAs present) for use when only "RSerr" is chosen; weights and offsets should not be used in the \code{lm} object}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw},
expected to be row-standardised (W-style)}
  \item{zero.policy}{default \code{attr(listw, "zero.policy")} as set when \code{listw} was created, if attribute not set, use global option value; if TRUE assign zero to the lagged value of zones without
neighbours, if FALSE assign NA}
  \item{test}{a character vector of tests requested chosen from RSerr, RSlag,
adjRSerr, adjRSlag, SARMA; test="all" computes all the tests.}
  \item{spChk}{should the data vector names be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use \code{get.spChkOption()}}
  \item{naSubset}{default TRUE to subset listw object for omitted observations in model object (this is a change from earlier behaviour, when the \code{model$na.action} component was ignored, and the listw object had to be subsetted by hand)}
  \item{x, object}{object to be printed}
  \item{p.adjust.method}{a character string specifying the probability value adjustment (see \code{\link{p.adjust}}) for multiple tests, default "none"}
  \item{digits}{minimum number of significant digits to be used for most numbers}
  \item{...}{printing arguments to be passed through}
}

\details{
The two types of dependence are for spatial lag \eqn{\rho}{rho} and spatial error \eqn{\lambda}{lambda}:

\deqn{
\mathbf{y} = \mathbf{X \beta} + \rho \mathbf{W_{(1)} y} + \mathbf{u},
}{y = X beta + rho W1 y + u}
\deqn{
\mathbf{u} = \lambda \mathbf{W_{(2)} u} + \mathbf{e}
}{u = lambda W2 u + e}

where \eqn{\mathbf{e}}{e} is a well-behaved, uncorrelated error
term. Tests for a missing spatially lagged dependent variable test
that \eqn{\rho = 0}{rho = 0}, tests for spatial autocorrelation of
the error \eqn{\mathbf{u}}{u} test whether \eqn{\lambda = 0}{lambda =
0}. \eqn{\mathbf{W}}{W} is a spatial weights matrix; for the tests used
here they are identical.

If factors or ordered factors are among explanatory variables in the regression, the coding of the contrasts affects the outcomes of significance tests.
}

\value{
A list of class \code{RStestlist} of \code{htest} objects, each with:
  \item{statistic}{the value of the Rao's score (a.k.a Lagrange multiplier) test.}
  \item{parameter}{number of degrees of freedom}
  \item{p.value}{the p-value of the test.}
  \item{method}{a character string giving the method used.}
  \item{data.name}{a character string giving the name(s) of the data.}
}
\references{Anselin, L. 1988 Spatial econometrics: methods and
models. (Dordrecht: Kluwer); Anselin, L., Bera, A. K., Florax, R. and
Yoon, M. J. 1996 Simple diagnostic tests for spatial dependence. Regional
Science and Urban Economics, 26, 77--104 \doi{10.1016/0166-0462(95)02111-6}; 
Malabika Koley (2024) Specification Testing under General Nesting Spatial 
Model, \url{https://sites.google.com/view/malabikakoley/research}.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no} and Andrew Bernat}

\seealso{\code{\link{lm}}, \code{\link{SD.RStests}}}

\examples{
data(oldcol)
oldcrime.lm <- lm(CRIME ~ HOVAL + INC, data = COL.OLD)
summary(oldcrime.lm)
lw <- nb2listw(COL.nb)
res <- lm.RStests(oldcrime.lm, listw=lw, test="all")
summary(res)
if (require("spatialreg", quietly=TRUE)) {
  oldcrime.slx <- lmSLX(CRIME ~ HOVAL + INC, data = COL.OLD, listw=lw)
  summary(lm.RStests(oldcrime.slx, listw=lw, test=c("adjRSerr", "adjRSlag")))
}
run <- require("codingMatrices", quietly=TRUE)
if (run) {
COL.OLD$fEW <- factor(COL.OLD$EW)
COL.OLD$fDISCBD <- ordered(cut(COL.OLD$DISCBD, c(0, 1.5, 3, 4.5, 6)))
f <- formula(CRIME ~ INC + HOVAL + fDISCBD*fEW)
lw <- nb2listw(COL.nb, style="W")
# default codings
summary(lm.RStests(lm(f, data=COL.OLD, contrasts=list(fDISCBD="contr.poly", fEW="contr.treatment")),
 lw, test=c("adjRSerr", "adjRSlag")))
}
if (run) {
# use codingMatrices::code_diff for ordered factor
summary(lm.RStests(lm(f, data=COL.OLD, contrasts=list(fDISCBD="code_diff", fEW="contr.treatment")),
 lw, test=c("adjRSerr", "adjRSlag")))
}
if (run) {
# use codingMatrices::code_control for both
summary(lm.RStests(lm(f, data=COL.OLD, contrasts=list(fDISCBD="code_control", fEW="code_control")),
 lw, test=c("adjRSerr", "adjRSlag")))
}
}
\keyword{spatial}
