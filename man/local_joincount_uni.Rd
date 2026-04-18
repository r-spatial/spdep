\name{local_joincount_uni}
\alias{local_joincount_uni}
\title{Calculate the local univariate join count}
\usage{
local_joincount_uni(
  fx,
  chosen,
  listw,
  alternative = "two.sided",
  nsim = 199,
  iseed = NULL,
  ties.method = "average",
  no_repeat_in_row=FALSE
)
}
\arguments{
\item{fx}{a factor with two levels; use of an ordered factor is not well understood.}

\item{chosen}{a scalar character containing the level of \code{fx} that should be considered the observed value (1).}

\item{listw}{a listw object containing binary weights created, for example, with \code{nbwlistw(nb, style = "B")}}

\item{alternative}{default \code{"greater"}. One of \code{"less"} or \code{"greater"}.}

\item{nsim}{the number of conditional permutation simulations}

\item{iseed}{default NULL, used to set the seed; the output will only be reproducible if the count of CPU cores across which computation is distributed is the same}

\item{ties.method}{default \code{"average"} passed through to \code{rank}, can take values accepted by \code{rank}: \code{c("average", "first", "last", "random", "max", "min")}, see \code{\link{rank}}}

\item{no_repeat_in_row}{default \code{FALSE}, if \code{TRUE}, sample conditionally in each row without replacements to avoid duplicate values, \url{https://github.com/r-spatial/spdep/issues/124}}

}
\value{
a \code{data.frame}  with class \dQuote{local-jc-uni} with six columns \code{BB} (observed BB - neighbour same as focus), \code{Pr()} (pseudo-p from punif rank), \code{sim_rank} (simulation rank with current ties.method), \code{p_sim_pysal_ge} (pseudo-p fromPySAL esda, greater than or equal to observed BB), \code{p_sim_pysal_gt} (pseudo-p fromPySAL esda, greater than observed BB), \code{largereq} (count of simulated values greater than or equal to observed BB) after folding, \code{olarger} (count of simulated values greater than observed BB) before folding, \code{olargereq} (count of simulated values greater than or equal to observed BB) before folding,and number of rows equal to the length of \code{x}.
}
\description{
The univariate local join count statistic is used to identify clusters of rarely occurring binary variables. The binary variable of interest should occur less
than half of the time.
}
\details{
The local join count statistic requires a binary weights list which can be generated with \code{nb2listw(nb, style = "B")}. Additionally, ensure that the binary variable of interest is rarely occurring in no more than half of observations.

P-values are estimated using a conditional permutation approach. This creates a reference distribution from which the observed statistic is compared. For more see \href{https://raw.githubusercontent.com/GeoDaCenter/GeoDaCenter.github.io/refs/heads/update-1.22/glossary.html#ppvalue}{Geoda Glossary}.

The pseudo-p-values returned by freestanding Geoda and PySAL esda correspond to \code{"res_min$p_sim_pysal_ge"} and ranked pseudo-p-values with \code{tied.method="min"}  (equivalently \code{"last"}; \code{"res_min$p_sim_pysal_gt"} corresponds to ranked pseudo-p-values with \code{tied.method="max"} (equivalently \code{"first"}).
}
\examples{
data(oldcol)
fx <- as.factor(ifelse(COL.OLD$CRIME < 35, "low-crime", "high-crime"))
listw <- nb2listw(COL.nb, style = "B")
set.seed(1)
res_min <- local_joincount_uni(fx, chosen = "high-crime", listw, nsim=999,
 alternative="two.sided", ties.method="min")
cor(res_min[,2], res_min$p_sim_pysal_ge, use="complete.obs")
res_max <- local_joincount_uni(fx, chosen = "high-crime", listw, nsim=999,
 alternative="two.sided", ties.method="max")
cor(res_max[,2], res_max$p_sim_pysal_gt, use="complete.obs")
}

\references{
Anselin, L., & Li, X. (2019). Operational Local Join Count Statistics for Cluster Detection. Journal of geographical systems, 21(2), 189–210. \doi{10.1007/s10109-019-00299-x}
}
\author{Josiah Parry \email{josiah.parry@gmail.com}}
