% Copyright 2003 by Roger S. Bivand
\name{anova.sarlm}
\alias{anova.sarlm}
\title{Comparison of simultaneous autoregressive models}
\description{
  One of a number of tools for comparing simultaneous autoregressive models, in particular nested models. The function is based on \code{anova.lme()} for comparing linear mixed models, and follows that function in using the "anova" generic name.
}
\usage{
\method{anova}{sarlm}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{\code{object} is of class \code{sarlm}}
  \item{\dots}{other objects of class \code{sarlm} or class \code{lm}}
}
\details{
  If successive models have different numbers of degrees of freedom, a likelihood ratio test will be performed between them. It is important to recall that tests apply to nested models, and this function at least attempts to make sure that the response variable in the models being compared has the same name. Useless results can still be generated when incomparable models are compared, it being the responsibility of the user to check.
}
\value{
  The function returns a data frame printed by default functions
}

\author{Roger Bivand \email{Roger.Bivand@nhh.no}}
 
\seealso{\code{\link{LR.sarlm}}, \code{\link{AIC}}}

\examples{
if (require(rgdal, quietly=TRUE)) {
example(columbus, package="spData")
lm.mod <- lm(CRIME ~ HOVAL + INC, data=columbus)
lag <- lagsarlm(CRIME ~ HOVAL + INC, data=columbus, nb2listw(col.gal.nb))
mixed <- lagsarlm(CRIME ~ HOVAL + INC, data=columbus, nb2listw(col.gal.nb),
  type="mixed")
error <- errorsarlm(CRIME ~ HOVAL + INC, data=columbus, nb2listw(col.gal.nb))
LR.sarlm(mixed, error)
anova(lag, lm.mod)
anova(lag, error, mixed)
AIC(lag, error, mixed)
}
}
\keyword{spatial}
