% Copyright 2014 by Roger S. Bivand, Virgilio Gómez-Rubio
\encoding{latin1}
\name{lee}
\alias{lee}
\title{Compute Lee's statistic}
\description{
  A simple function to compute Lee's L statistic for bivariate spatial data;
%, called by \code{moran.test} and \code{moran.mc};
\deqn{L(x,y) = \frac{n}{\sum_{i=1}^{n}(\sum_{j=1}^{n}w_{ij})^2}
\frac{\sum_{i=1}^{n}(\sum_{j=1}^{n}w_{ij}(x_i-\bar{x})) ((\sum_{j=1}^{n}w_{ij}(y_j-\bar{y}))}{\sqrt{\sum_{i=1}^{n}(x_i - \bar{x})^2} \sqrt{\sum_{i=1}^{n}(y_i - \bar{y})^2}}
}{L(x,y) = (n sum_i (sum_j w_ij (x_i - xbar)) (sum_j w_ij (y_j - ybar))) / (S2 sqrt(sum_i (x_i - xbar)^2)) sqrt(sum_i (x_i - xbar)^2))}
}
\usage{
lee(x, y, listw, n, S2, zero.policy=NULL, NAOK=FALSE)
}
\arguments{
  \item{x}{a numeric vector the same length as the neighbours list in listw}
  \item{y}{a numeric vector the same length as the neighbours list in listw}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{n}{number of zones}
  \item{S2}{Sum of squared sum of weights by rows.}
  \item{zero.policy}{default NULL, use global option value; if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA}
  \item{NAOK}{if 'TRUE' then any 'NA' or 'NaN' or 'Inf' values in x are passed on to the foreign function.  If 'FALSE', the presence of 'NA' or 'NaN' or 'Inf' values is regarded as an error.}
}

\value{
a list of
  \item{L}{Lee's L statistic}
  \item{local L}{Lee's local L statistic}
}
\references{Lee (2001). Developing a bivariate spatial association measure:
An integration of Pearson's r and Moran's I. J Geograph Syst  3: 369-385}
\author{Roger Bivand and Virgiio Gómez-Rubio \email{Virgilio.Gomez@uclm.es}}

\seealso{\code{\link{lee.mc}}}%, \code{\link{lee.mc}}}

\examples{
data(boston, package="spData")
lw<-nb2listw(boston.soi)

x<-boston.c$CMEDV
y<-boston.c$CRIM
z<-boston.c$RAD

Lxy<-lee(x, y, lw, length(x), zero.policy=TRUE)
Lxz<-lee(x, z, lw, length(x), zero.policy=TRUE)

}
\keyword{spatial}
