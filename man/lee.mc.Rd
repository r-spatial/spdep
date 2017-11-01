% Copyright 2014 by Roger S. Bivand, Virgilio Gómez-Rubio
\encoding{latin1}
\name{lee.mc}
\alias{lee.mc}
\title{Permutation test for Lee's L statistic}
\description{
 A permutation test for Lee's L  statistic calculated by using nsim random permutations of x and y for the given spatial weighting scheme, to establish the rank of the observed statistic in relation to the nsim simulated values. 
}
\usage{
lee.mc(x, y, listw, nsim, zero.policy=NULL, alternative="greater",
 na.action=na.fail, spChk=NULL, return_boot=FALSE)
}
\arguments{
  \item{x}{a numeric vector the same length as the neighbours list in listw}
  \item{y}{a numeric vector the same length as the neighbours list in listw}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{nsim}{number of permutations}
  \item{zero.policy}{default NULL, use global option value; if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA}
  \item{alternative}{a character string specifying the alternative hypothesis, must be one of "greater" (default), or "less".}
  \item{na.action}{a function (default \code{na.fail}), can also be \code{na.omit} or \code{na.exclude} - in these cases the weights list will be subsetted to remove NAs in the data. It may be necessary to set zero.policy to TRUE because this subsetting may create no-neighbour observations. Note that only weights lists created without using the glist argument to \code{nb2listw} may be subsetted. \code{na.pass} is not permitted because it is meaningless in a permutation test.}
  \item{spChk}{should the data vector names be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use \code{get.spChkOption()}}
  \item{return_boot}{return an object of class \code{boot} from the equivalent permutation bootstrap rather than an object of class \code{htest}}
}

\value{
A list with class \code{htest} and \code{mc.sim} containing the following components:
  \item{statistic}{the value of the observed Lee's L.}
  \item{parameter}{the rank of the observed Lee's L.}
  \item{p.value}{the pseudo p-value of the test.}
  \item{alternative}{a character string describing the alternative hypothesis.}
  \item{method}{a character string giving the method used.}
  \item{data.name}{a character string giving the name(s) of the data, and the number of simulations.}
  \item{res}{nsim simulated values of statistic, final value is observed statistic}
}
\references{Lee (2001). Developing a bivariate spatial association measure:
An integration of Pearson's r and Moran's I. J Geograph Syst  3: 369-385}

\author{Roger Bivand, Virgilio Gómez-Rubio \email{Virgilio.Gomez@uclm.es} }

\seealso{\code{\link{lee}}}%, \code{\link{moran.test}}}

\examples{
data(boston, package="spData")
lw<-nb2listw(boston.soi)

x<-boston.c$CMEDV
y<-boston.c$CRIM

lee.mc(x, y, nsim=99, lw, zero.policy=TRUE, alternative="less")

#Test with missing values
x[1:5]<-NA
y[3:7]<-NA

lee.mc(x, y, nsim=99, lw, zero.policy=TRUE, alternative="less", 
   na.action=na.omit)

}
\keyword{spatial}
