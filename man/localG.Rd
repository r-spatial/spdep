% Copyright 2001 by Roger S. Bivand
\name{localG}
\alias{localG}
\title{G and Gstar local spatial statistics}
\description{
The local spatial statistic G is calculated for each zone based on the
spatial weights object used. The value returned is a Z-value, and may be
used as a diagnostic tool. High positive values indicate the posibility
of a local cluster of high values of the variable being analysed, very
low relative values a similar cluster of low values. For inference,
a Bonferroni-type test is suggested in the references, where tables of
critical values may be found (see also details below).
}
\usage{
localG(x, listw, zero.policy=NULL, spChk=NULL, return_internals=FALSE)
}
\arguments{
  \item{x}{a numeric vector the same length as the neighbours list in listw}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{zero.policy}{default NULL, use global option value; if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA}
  \item{spChk}{should the data vector names be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use \code{get.spChkOption()}}
  \item{return_internals}{default FALSE, if TRUE, return internal values of G, EI and VG as as attribute matrix}
}
\details{
If the neighbours member of listw has a "self.included" attribute set
to TRUE, the Gstar variant, including the self-weight \eqn{w_{ii} > 0},
is calculated and returned.  The returned vector will have a "gstari"
attribute set to TRUE.  Self-weights can be included by using the
\code{include.self} function before converting
the neighbour list to a spatial weights list with \code{nb2listw} as
shown below in the example.

The critical values of the statistic under assumptions given in the
references for the 95th percentile are for n=1: 1.645, n=50: 3.083,
n=100: 3.289, n=1000: 3.886.
}
\value{
 A vector of G or Gstar values, with attributes "gstari" set to TRUE or
FALSE, "call" set to the function call, and class "localG".
}
\references{Ord, J. K. and Getis, A. 1995 Local spatial autocorrelation
statistics: distributional issues and an application. \emph{Geographical
Analysis}, 27, 286--306; Getis, A. and Ord, J. K. 1996 Local spatial
statistics: an overview. In P. Longley and M. Batty (eds) \emph{Spatial
analysis: modelling in a GIS environment} (Cambridge: Geoinformation
International), 261--277.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

%\seealso{\code{\link{localmoran}}}

\examples{
data(getisord, package="spData")
xycoords <- cbind(xyz$x, xyz$y)
nb30 <- dnearneigh(xycoords, 0, 30)
G30 <- localG(xyz$val, nb2listw(nb30, style="B"))
G30[length(xyz$val)-136]
nb60 <- dnearneigh(xycoords, 0, 60)
G60 <- localG(xyz$val, nb2listw(nb60, style="B"))
G60[length(xyz$val)-136]
nb90 <- dnearneigh(xycoords, 0, 90)
G90 <- localG(xyz$val, nb2listw(nb90, style="B"))
G90[length(xyz$val)-136]
nb120 <- dnearneigh(xycoords, 0, 120)
G120 <- localG(xyz$val, nb2listw(nb120, style="B"))
G120[length(xyz$val)-136]
nb150 <- dnearneigh(xycoords, 0, 150)
G150 <- localG(xyz$val, nb2listw(nb150, style="B"))
G150[length(xyz$val)-136]
brks <- seq(-5,5,1)
cm.col <- cm.colors(length(brks)-1)
image(x, y, t(matrix(G30, nrow=16, ncol=16, byrow=TRUE)),
  breaks=brks, col=cm.col, asp=1)
text(xyz$x, xyz$y, round(G30, digits=1), cex=0.7)
polygon(c(195,225,225,195), c(195,195,225,225), lwd=2)
title(main=expression(paste("Values of the ", G[i], " statistic")))
G30s <- localG(xyz$val, nb2listw(include.self(nb30),
 style="B"))
cat("value according to Getis and Ord's eq. 14.2, p. 263 (1996)\n")
G30s[length(xyz$val)-136]
cat(paste("value given by Getis and Ord (1996), p. 267",
  "(division by n-1 rather than n \n in variance)\n"))
G30s[length(xyz$val)-136] *
  (sqrt(sum(scale(xyz$val, scale=FALSE)^2)/length(xyz$val)) /
  sqrt(var(xyz$val)))
image(x, y, t(matrix(G30s, nrow=16, ncol=16, byrow=TRUE)),
  breaks=brks, col=cm.col, asp=1)
text(xyz$x, xyz$y, round(G30s, digits=1), cex=0.7)
polygon(c(195,225,225,195), c(195,195,225,225), lwd=2)
title(main=expression(paste("Values of the ", G[i]^"*", " statistic")))
}
\keyword{spatial}
