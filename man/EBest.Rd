% Copyright 2003-6 by Roger S. Bivand
\name{EBest}
\alias{EBest}
\title{Global Empirical Bayes estimator}
\description{
  The function computes global empirical Bayes estimates for rates "shrunk" to the overall mean.
}
\usage{
EBest(n, x, family="poisson")
}
\arguments{
  \item{n}{a numeric vector of counts of cases}
  \item{x}{a numeric vector of populations at risk}
  \item{family}{either "poisson" for rare conditions or "binomial" for non-rare conditions}
}
\details{
  Details of the implementation for the "poisson" family are to be found in Marshall, p. 284--5, and Bailey and Gatrell p. 303--306 and exercise 8.2, pp. 328--330. For the "binomial" family, see Martuzzi and Elliott (implementation by Olaf Berke).
}
\value{
A data frame with two columns:
  \item{raw}{a numerical vector of raw (crude) rates}
  \item{estmm}{a numerical vector of empirical Bayes estimates}

and a \code{parameters} attribute list with components:

  \item{a}{global method of moments phi value}
  \item{m}{global method of moments gamma value}
}
\references{Marshall R M (1991) Mapping disease and mortality rates using Empirical Bayes Estimators, Applied Statistics, 40, 283--294; Bailey T, Gatrell A (1995) Interactive Spatial Data Analysis, Harlow: Longman, pp. 303--306, Martuzzi M, Elliott P (1996) Empirical Bayes estimation of small area prevalence of non-rare conditions, Statistics in Medicine 15, 1867--1873.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no} and Olaf Berke, Population Medicine, OVC, University of Guelph, CANADA}

\seealso{\code{\link{EBlocal}}, \code{\link{probmap}}, \code{\link{EBImoran.mc}}}

\examples{
if (require(rgdal, quietly=TRUE)) {
example(auckland, package="spData")
res <- EBest(auckland$M77_85, 9*auckland$Und5_81)
attr(res, "parameters")
if (require(classInt, quietly=TRUE)) {
cols <- grey(6:2/7)
cI <- classIntervals(res$estmm*1000, style="fixed",
 fixedBreaks=c(-Inf,2,2.5,3,3.5,Inf))
fcI <- findColours(cI, pal=grey(6:2/7))
plot(auckland, col=fcI)
legend("bottomleft", fill=attr(fcI, "palette"),
 legend=names(attr(fcI, "table")), bty="n")
title(main="Global moment estimator of infant mortality per 1000 per year")
}
}
data(huddersfield, package="spData")
res <- EBest(huddersfield$cases, huddersfield$total, family="binomial")
round(res[,1:2],4)*100
}
\keyword{spatial}
