% Copyright 2002 by Roger S. Bivand
\name{sp.correlogram}
\alias{sp.correlogram}
\alias{plot.spcor}
\alias{print.spcor}
\title{Spatial correlogram}
\description{
  Spatial correlograms for Moran's I and the autocorrelation coefficient, with
print and plot helper functions.
}
\usage{
sp.correlogram(neighbours, var, order = 1, method = "corr",
 style = "W", randomisation = TRUE, zero.policy = NULL, spChk=NULL)
\method{plot}{spcor}(x, main, ylab, ylim, ...)
\method{print}{spcor}(x, p.adj.method="none", ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{neighbours}{an object of class \code{nb}}
  \item{var}{a numeric vector}
  \item{order}{maximum lag order}
  \item{method}{"corr" for correlation, "I" for Moran's I, "C" for Geary's C}
  \item{style}{\code{style} can take values W, B, C, and S}
  \item{randomisation}{variance of I or C calculated under the assumption of randomisation, if FALSE normality}
  \item{zero.policy}{default NULL, use global option value; if FALSE stop with error for any empty neighbour sets, 
if TRUE permit the weights list to be formed with zero-length weights vectors}
  \item{spChk}{should the data vector names be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use \code{get.spChkOption()}}
  \item{x}{an object from \code{sp.correlogram()} of class \code{spcor}}
  \item{p.adj.method}{correction method as in \code{p.adjust}}
  \item{main}{an overall title for the plot}
  \item{ylab}{a title for the y axis}
  \item{ylim}{the y limits of the plot}
  \item{\dots}{further arguments passed through}
}
\details{
  The print function also calculates the standard deviates of Moran's I or Geary's C and a two-sided probability value, optionally using \code{p.adjust} to correct by the nymber of lags. The plot function plots a bar from the estimated Moran's I, or Geary's C value to +/- twice the square root of its variance (in previous releases only once, not twice). The table includes the count of included observations in brackets after the lag order. Care needs to be shown when interpreting results for few remaining included observations as lag order increases.
}

\value{
  returns a list of class \code{spcor}:
  \item{res}{for "corr" a vector of values; for "I", a matrix of estimates of "I", expectations, and variances}
  \item{method}{"I" or "corr"}
  \item{cardnos}{list of tables of neighbour cardinalities for the lag orders used}
  \item{var}{variable name}
}
\references{Cliff, A. D., Ord, J. K. 1981 \emph{Spatial processes}, Pion, pp. 
118--122, Martin, R. L., Oeppen, J. E. 1975 The identification of regional
forecasting models using space-time correlation functions, \emph{Transactions
of the Institute of British Geographers}, 66, 95--118.}
\author{Roger Bivand, \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{nblag}}, \code{\link{moran}}, \code{\link[stats]{p.adjust}}}

\examples{
if (require(rgdal, quietly=TRUE)) {
example(nc.sids, package="spData")
ft.SID74 <- sqrt(1000)*(sqrt(nc.sids$SID74/nc.sids$BIR74) +
  sqrt((nc.sids$SID74+1)/nc.sids$BIR74))
tr.SIDS74 <- ft.SID74*sqrt(nc.sids$BIR74)
cspc <- sp.correlogram(ncCC89_nb, tr.SIDS74, order=8, method="corr",
 zero.policy=TRUE)
print(cspc)
plot(cspc)
Ispc <- sp.correlogram(ncCC89_nb, tr.SIDS74, order=8, method="I",
 zero.policy=TRUE)
print(Ispc)
print(Ispc, "bonferroni")
plot(Ispc)
Cspc <- sp.correlogram(ncCC89_nb, tr.SIDS74, order=8, method="C",
 zero.policy=TRUE)
print(Cspc)
print(Cspc, "bonferroni")
plot(Cspc)
drop.no.neighs <- !(1:length(ncCC89_nb) \%in\% which(card(ncCC89_nb) == 0))
sub.ncCC89.nb <- subset(ncCC89_nb, drop.no.neighs)
plot(sp.correlogram(sub.ncCC89.nb, subset(tr.SIDS74,  drop.no.neighs),
 order=8, method="corr"))
}
}
\keyword{spatial}
