% Copyright 2018-01 by Rene Westerholt
\name{LOSH}
\alias{LOSH}
\title{Local spatial heteroscedasticity}
  \description{
Local spatial heteroscedasticity is calculated for each location based on the spatial weights object used. The statistic is:
    \deqn{H_i = \frac{\sum_j^n w_{ij} \cdot |e_j|^a}{h_1 \cdot \sum_j^n w_{ij}}} with \deqn{e_j = x_j - \bar{x}_j} and \deqn{\bar{x}_j = \frac{\sum_k^n w_{jk} \cdot x_k}{\sum_k^n w_{jk}}}
    Its expectation and variance are given in Ord & Getis (2012). The exponent \emph{a} allows for investigating different types of mean dispersal.
  }
\usage{
  LOSH(x, listw, a=2, var_hi=TRUE, zero.policy=attr(listw, "zero.policy"),
 na.action=na.fail, spChk=NULL)
}
\arguments{
  \item{x}{a numeric vector of the same length as the neighbours list in listw}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{a}{the exponent applied to the local residuals; the default value of 2 leads to a measure of heterogeneity in the spatial variance}
  \item{var_hi}{default TRUE, the moments and the test statistics are calculated for each location; if FALSE, only the plain LOSH measures, \eqn{\bar{x}_i} and \eqn{e_i} are calculated}
  \item{zero.policy}{default \code{attr(listw, "zero.policy")} as set when \code{listw} was created, if attribute not set, use global option value; if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA}
  \item{na.action}{a function (default \code{na.fail}), can also be \code{na.omit} or \code{na.exclude} - in these cases the weights list will be subsetted to remove NAs in the data. It may be necessary to set zero.policy to TRUE because this subsetting may create no-neighbour observations. Note that only weights lists created without using the glist argument to \code{nb2listw} may be subsetted. If \code{na.pass} is used, zero is substituted for NA values in calculating the spatial lag. (Note that na.exclude will only work properly starting from R 1.9.0, na.omit and na.exclude assign the wrong classes in 1.8.*)}
    \item{spChk}{should the data vector names be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use \code{get.spChkOption()}}
  }
    \details{
    In addition to the LOSH measure, the values returned include local spatially weighted mean values \eqn{\bar{x}_i} and local residuals \eqn{e_i} estimated about these means. These values facilitate the interpretation of LOSH values. Further, if specified through \code{var_hi}, the statistical moments and the test statistics as proposed by Ord & Getis (2012) are also calculated and returned.
    }
  
  \value{
    \item{Hi}{LOSH statistic}
    \item{E.Hi}{(optional) expectation of LOSH}
    \item{Var.Hi}{(optional) variance of LOSH}
    \item{Z.Hi}{(optional) the approximately Chi-square distributed test statistics}
    \item{x_bar_i}{local spatially weighted mean values}
    \item{ei}{residuals about local spatially weighted mean values}
  }
  
  \references{Ord, J. K., & Getis, A. 2012. Local spatial heteroscedasticity (LOSH),
    The Annals of Regional Science, 48 (2), 529--539.}
  \author{René Westerholt \email{rene.westerholt@tu-dortmund.de}}
  
  \seealso{\code{\link{LOSH.cs}}, \code{\link{LOSH.mc}}}
  
  \examples{
    data(boston, package="spData")
    resLOSH <- LOSH(boston.c$NOX, nb2listw(boston.soi))
    hist(resLOSH[,"Hi"])
    mean(resLOSH[,"Hi"])
  }
  \keyword{spatial}
