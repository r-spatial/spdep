% Copyright 2002-9 by Roger S. Bivand
\name{LR.sarlm}
\alias{LR.sarlm}
\alias{LR1.sarlm}
\alias{Wald1.sarlm}
\alias{Hausman.test}
\alias{Hausman.test.sarlm}
\alias{Hausman.test.gmsar}
\alias{logLik.sarlm}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Likelihood ratio test}
\description{
  The \code{LR.sarlm()} function provides a likelihood ratio test for objects for which a \code{logLik()} function exists for their class, or for objects of class \code{logLik}. \code{LR1.sarlm()} and \code{Wald1.sarlm()} are used internally in \code{summary.sarlm()}, but may be accessed directly; they report the values respectively of LR and Wald tests for the absence of spatial dependence in spatial lag or error models. The spatial Hausman test is available for models fitted with \code{errorsarlm} and \code{GMerrorsar}.
}
\usage{
LR.sarlm(x, y)
\method{logLik}{sarlm}(object, ...)
LR1.sarlm(object)
Wald1.sarlm(object)
\method{Hausman.test}{sarlm}(object, ..., tol=NULL)
\method{Hausman.test}{gmsar}(object, ..., tol=NULL)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{a \code{logLik} object or an object for which a \code{logLik()} function exists}
  \item{y}{a \code{logLik} object or an object for which a \code{logLik()} function exists}
  \item{object}{a \code{sarlm} object from \code{lagsarlm()} or \code{errorsarlm()}}
  \item{\dots}{further arguments passed to or from other methods}
  \item{tol}{\code{tol} argument passed to \code{solve}, default NULL}
}
\value{
  The tests return objects of class \code{htest} with:
  \item{statistic}{value of statistic}
  \item{parameter}{degrees of freedom}
  \item{p.value}{Probability value}
  \item{estimate}{varies with test}
  \item{method}{description of test method}

  \code{logLik.sarlm()} returns an object of class \code{logLik}
  \code{LR1.sarlm}, \code{Hausman.sarlm} and \code{Wald1.sarlm} returm objects of class \code{htest}
}

\note{The numbers of degrees of freedom returned by \code{logLik.sarlm()} include nuisance parameters, that is the number of regression coefficients, plus sigma, plus spatial parameter esitmate(s).}

\references{LeSage J and RK Pace (2009) Introduction to Spatial Econometrics. CRC Press, Boca Raton, pp. 61--63; Pace RK and LeSage J (2008) A spatial Hausman test. \emph{Economics Letters} 101, 282--284.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{logLik.lm}}, \code{\link{anova.sarlm}}}

\examples{
if (require(rgdal, quietly=TRUE)) {
example(columbus, package="spData")
mixed <- lagsarlm(CRIME ~ HOVAL + INC, data=columbus, nb2listw(col.gal.nb),
  type="mixed")
error <- errorsarlm(CRIME ~ HOVAL + INC, data=columbus, nb2listw(col.gal.nb))
LR.sarlm(mixed, error)
Hausman.test(error)
}
}
\keyword{spatial}
