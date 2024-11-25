\name{localC}
\alias{localC}
\alias{localC.default}
\alias{localC.formula}
\alias{localC.list}
\alias{localC.matrix}
\alias{localC.data.frame}
\alias{localC_perm}
\alias{localC_perm.default}
\alias{localC_perm.formula}
\title{Compute Local Geary statistic}
\usage{
localC(x, ..., zero.policy=NULL)

\method{localC}{default}(x, listw, ..., zero.policy=attr(listw, "zero.policy"))

\method{localC}{formula}(formula, data, listw, ..., zero.policy=attr(listw, "zero.policy"))

\method{localC}{list}(x, listw, ..., zero.policy=attr(listw, "zero.policy"))

\method{localC}{matrix}(x, listw, ..., zero.policy=attr(listw, "zero.policy"))

\method{localC}{data.frame}(x, listw, ..., zero.policy=attr(listw, "zero.policy"))

localC_perm(x, ..., zero.policy=NULL, iseed=NULL, no_repeat_in_row=FALSE)

\method{localC_perm}{default}(x, listw, nsim = 499, alternative = "two.sided", ...,
 zero.policy=attr(listw, "zero.policy"), iseed=NULL, no_repeat_in_row=FALSE)

\method{localC_perm}{formula}(formula, data, listw, nsim = 499,
 alternative = "two.sided", ..., zero.policy=attr(listw, "zero.policy"), iseed=NULL,
 no_repeat_in_row=FALSE)
}
\arguments{
\item{x}{a numeric vector, numeric matrix, or list. See details for more.}

\item{formula}{A one-sided formula determining which variables to be used.}

\item{listw}{a \code{listw} object created for example by \code{nb2listw}.}

\item{data}{Used when a formula is provided. A matrix or data frame containing the variables in the formula \code{formula}.}

\item{nsim}{The number of simulations to be used for permutation test.}

\item{alternative}{A character defining the alternative hypothesis. Must be one of \code{"two.sided"}, \code{"less"} or \code{"greater"}.}
\item{...}{other arguments passed to methods.}
\item{zero.policy}{default \code{attr(listw, "zero.policy")} as set when \code{listw} was created, if attribute not set, use global option value; if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA.}
\item{iseed}{default NULL, used to set the seed;the output will only be reproducible if the count of CPU cores across which computation is distributed is the same}
\item{no_repeat_in_row}{default \code{FALSE}, if \code{TRUE}, sample conditionally in each row without replacements to avoid duplicate values, \url{https://github.com/r-spatial/spdep/issues/124}}
}
\description{
The Local Geary is a local adaptation of Geary's C statistic of spatial autocorrelation. The Local Geary uses squared differences to measure dissimilarity unlike the Local Moran. Low values of the Local Geary indicate positive spatial autocorrelation and large refers to negative spatial autocorrelation.

Inference for the Local Geary is based on a permutation approach which compares the observed value to the reference distribution under spatial randomness. \code{localC_perm()} returns a pseudo p-value. This is not an analytical p-value and is based on the number of permutations and as such should be used with care.
}
\details{
The Local Geary can be extended to a multivariate context. When \code{x} is a numeric vector, the univariate Local Geary will be calculated. To calculate the multivariate Local Moran provide either a list or a matrix. When \code{x} is a list, each element must be a numeric vector of the same length and of the same length as the neighbours in \code{listw}. In the case that \code{x} is a matrix the number of rows must be the same as the length of the neighbours in \code{listw}.

While not required in the univariate context, the standardized Local Geary is calculated. The multivariate Local Geary is \emph{always} standardized.

The univariate Local Geary is calculated as \eqn{c_i = \sum_j w_{ij}(x_i - x_j)^2} and the multivariate Local Geary is calculated as \eqn{c_{k,i} = \sum_{v=1}^{k} c_{v,i}} as described in Anselin (2019).
}

\value{
A numeric vector containing Local Geary statistic with attribute \code{pseudo-p} when \code{localC_perm()} is used. \code{pseudo-p} is an 8 column matrix containing

  \item{E.Ci}{expectation of the Local Geary statistic based on permutation sample}
  \item{Var.Ci}{variance of Local Geary based on permutation sample}
  \item{Z.Ci}{standard deviate of Local Geary based on permutation sample}
  \item{Pr()}{p-value of Local Geary statistic using \code{pnorm()} using standard deviates based on permutation sample means and standard deviations}
  \item{Pr() Sim}{\code{rank()} and \code{punif()} of observed statistic rank for [0, 1] p-values using \code{alternative=}}
  \item{Pr(folded) Sim}{the simulation folded [0, 0.5] range ranked p-value (based on \url{https://github.com/pysal/esda/blob/4a63e0b5df1e754b17b5f1205b8cadcbecc5e061/esda/crand.py#L211-L213})}
  \item{Skewness}{the output of \code{e1071::skewness()} for the permutation samples underlying the standard deviates}
  \item{Kurtosis}{the output of \code{e1071::kurtosis()} for the permutation samples underlying the standard deviates}
}

\examples{
orig <- spData::africa.rook.nb
listw <- nb2listw(orig)
x <- spData::afcon$totcon

(A <- localC(x, listw))
listw1 <- nb2listw(droplinks(sym.attr.nb(orig), 3, sym=TRUE), zero.policy=TRUE)
(A1 <- localC(x, listw1, zero.policy=FALSE))
(A2 <- localC(x, listw1, zero.policy=TRUE))


\dontrun{
library(reticulate)
use_python("/usr/bin/python", required = TRUE)
gp <- import("geopandas")
ps <- import("libpysal")
W <- listw2mat(listw)
w <- ps$weights$full2W(W, rownames(W))
w$transform <- "R"
esda <- import("esda")
lM <- esda$Moran_Local(x, w)
all.equal(unname(localmoran(x, listw, mlvar=FALSE)[,1]), c(lM$Is))
# confirm x and w the same
lC <- esda$Geary_Local(connectivity=w)$fit(scale(x))
# np$std missing ddof=1
n <- length(x)
D0 <- spdep:::geary.intern((x - mean(x)) / sqrt(var(x)*(n-1)/n), listw, n=n)
# lC components probably wrongly ordered https://github.com/pysal/esda/issues/192
o <- match(round(D0, 6), round(lC$localG, 6))
all.equal(c(lC$localG)[o], D0)
# simulation order not retained
lC$p_sim[o]
attr(C, "pseudo-p")[,6]
}
}
\references{
Anselin, L. (1995), Local Indicators of Spatial Association—LISA. Geographical Analysis, 27: 93-115. \doi{10.1111/j.1538-4632.1995.tb00338.x}

Anselin, L. (2019), A Local Indicator of Multivariate Spatial Association: Extending Gearys c. Geogr Anal, 51: 133-150. \doi{10.1111/gean.12164}
}
\author{
Josiah Parry \email{josiah.parry@gmail.com} and Roger Bivand
}
