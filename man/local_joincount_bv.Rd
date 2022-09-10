\name{local_joincount_bv}
\alias{local_joincount_bv}
\title{Calculate the local bivariate join count}
\usage{
local_joincount_bv(
  x,
  z,
  listw,
  nsim = 199,
  alternative = "two.sided"
)
}

\arguments{
\item{x}{a binary variable either numeric or logical}
\item{z}{a binary variable either numeric or logical with the same length as \code{x}}

\item{listw}{a listw object containing binary weights created, for example, with \code{nb2listw(nb, style = "B")}}

\item{nsim}{the number of conditional permutation simulations}

\item{alternative}{default \code{"greater"}. One of \code{"less"} or \code{"greater"}.}
}

\description{
The bivariate join count (BJC) evaluates event occurrences in predefined regions and tests if the co-occurrence of events deviates from complete spatial randomness.
}

\details{
  There are two cases that are evaluated in the bivariate join count. The first being in-situ colocation (CLC) where xi = 1 and zi = 1. The second is the general form of the bivariate join count (BJC) that is used when there is no in-situ colocation.

  The BJC case "is useful when x and z cannot occur in the same location, such as when x and z correspond to two different values of a single categorical variable" or "when x and z can co-locate, but do not" (Anselin and Li, 2019). Whereas the CLC case is useful in evaluating simultaneous occurrences of events.

  The local bivariate join count statistic requires a binary weights list which can be generated with \code{nb2listw(nb, style = "B")}.

  P-values are estimated using a conditional permutation approach. This creates a reference distribution from which the observed statistic is compared.
}

\value{
a \code{data.frame} with two columns \code{join_count} and \code{p_sim} and number of rows equal to the length of arguments \code{x}.
}


\references{
Anselin, L., & Li, X. (2019). Operational Local Join Count Statistics for Cluster Detection. Journal of geographical systems, 21(2), 189–210. \doi{10.1007/s10109-019-00299-x}
}

\author{Josiah Parry \email{josiah.parry@gmail.com}}
