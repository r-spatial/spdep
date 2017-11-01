% Copyright 2002-3 by Roger S. Bivand
\name{residuals.sarlm}
\alias{residuals.sarlm}
\alias{deviance.sarlm}
\alias{coef.sarlm}
\alias{vcov.sarlm}
\alias{fitted.sarlm}
\title{Access functions for spatial simultaneous autoregressive linear
model objects}
\description{
  Access functions for residuals, deviance, coefficients and fitted values 
from spatial simultaneous autoregressive linear model objects
}
\usage{
\method{residuals}{sarlm}(object, ...)
\method{deviance}{sarlm}(object, ...)
\method{coef}{sarlm}(object, ...)
\method{vcov}{sarlm}(object, ...)
\method{fitted}{sarlm}(object, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{\code{sarlm} object returned by \code{lagsarlm} or 
\code{errorsarlm}}
  \item{...}{further arguments passed through}
}

\value{
  Revelant vectors of numerical values.
}

\note{\code{fitted.sarlm()} returns the difference between 
\code{residuals()} for the same object and the response variable; 
\code{predict.sarlm()} returns a decomposition into trend and signal 
for the fit.}

\author{Roger Bivand, \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{errorsarlm}}, \code{\link{lagsarlm}}, 
\code{\link{predict.sarlm}}}

\keyword{spatial}
