% Copyright 2004-2010 by Roger S. Bivand
\name{choynowski}
\alias{choynowski}
\title{Choynowski probability map values}
\description{
  Calculates Choynowski probability map values.
}
\usage{
choynowski(n, x, row.names=NULL, tol = .Machine$double.eps^0.5, legacy=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{a numeric vector of counts of cases}
  \item{x}{a numeric vector of populations at risk}
  \item{row.names}{row names passed through to output data frame}
  \item{tol}{accumulate values for observed counts >= expected until value less than tol}
  \item{legacy}{default FALSE using vectorised alternating side \code{ppois} version, if true use original version written from sources and iterating down to \code{tol}}
}

\value{
  A data frame with columns:
  \item{pmap}{Poisson probability map values: probablility of getting a more ``extreme'' count than actually observed, one-tailed with less than expected and more than expected folded together}
  \item{type}{logical: TRUE if observed count less than expected}
}
\references{Choynowski, M (1959) Maps based on probabilities, Journal of
the American Statistical Association, 54, 385--388; Cressie, N, Read,
TRC (1985), Do sudden infant deaths come in clusters? Statistics and
Decisions, Supplement Issue 2, 333--349; Bailey T, Gatrell A (1995)
Interactive Spatial Data Analysis, Harlow: Longman, pp. 300--303.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{probmap}}}
\examples{
if (require(rgdal, quietly=TRUE)) {
example(auckland, package="spData")
res <- choynowski(auckland$M77_85, 9*auckland$Und5_81)
resl <- choynowski(auckland$M77_85, 9*auckland$Und5_81, legacy=TRUE)
all.equal(res, resl)
rt <- sum(auckland$M77_85)/sum(9*auckland$Und5_81)
ch_ppois_pmap <- numeric(length(auckland$Und5_81))
side <- c("greater", "less")
for (i in seq(along=ch_ppois_pmap)) {
  ch_ppois_pmap[i] <- poisson.test(auckland$M77_85[i], r=rt,
    T=(9*auckland$Und5_81[i]), alternative=side[(res$type[i]+1)])$p.value
}
all.equal(ch_ppois_pmap, res$pmap)

res1 <- probmap(auckland$M77_85, 9*auckland$Und5_81)
table(abs(res$pmap - res1$pmap) < 0.00001, res$type)
lt005 <- (res$pmap < 0.05) & (res$type)
ge005 <- (res$pmap < 0.05) & (!res$type)
cols <- rep("white", length(lt005))
cols[lt005] <- grey(2/7)
cols[ge005] <- grey(5/7)
plot(auckland, col=cols) 
legend("bottomleft", fill=grey(c(2,5)/7), legend=c("low", "high"), bty="n")
}
}
\keyword{spatial}
