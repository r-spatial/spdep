% Copyright 2014 by Roger S. Bivand,, Virgilio Gómez-Rubio
\encoding{latin1}
\name{lee.test}
\alias{lee.test}
\title{Lee's L  test for spatial autocorrelation}
\description{
Lee's L test for spatial autocorrelation using a spatial weights matrix in weights list form. The assumptions underlying the test are sensitive to the form of the graph of neighbour relationships and other factors, and results may be checked against those of \code{lee.mc} permutations. 
}
\usage{
lee.test(x, y, listw, zero.policy=NULL,
 alternative="greater", na.action=na.fail, spChk=NULL)
}
\arguments{
  \item{x}{a numeric vector the same length as the neighbours list in listw}
  \item{y}{a numeric vector the same length as the neighbours list in listw}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
%  \item{randomisation}{variance of I calculated under the assumption of randomisation, if FALSE normality}
  \item{zero.policy}{default NULL, use global option value; if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA}
  \item{alternative}{a character string specifying the alternative hypothesis, must be one of greater (default), less or two.sided.}
%  \item{rank}{logical value - default FALSE for continuous variables, if TRUE, uses the adaptation of Moran's I for ranks suggested by Cliff and Ord (1981, p. 46)}
  \item{na.action}{a function (default \code{na.fail}), can also be \code{na.omit} or \code{na.exclude} - in these cases the weights list will be subsetted to remove NAs in the data. It may be necessary to set zero.policy to TRUE because this subsetting may create no-neighbour observations. Note that only weights lists created without using the glist argument to \code{nb2listw} may be subsetted. If \code{na.pass} is used, zero is substituted for NA values in calculating the spatial lag}
  \item{spChk}{should the data vector names be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use \code{get.spChkOption()}}
%  \item{adjust.n}{default TRUE, if FALSE the number of observations is not adjusted for no-neighbour observations, if TRUE, the number of observations is adjusted}
}


\value{
A list with class \code{htest} containing the following components:
  \item{statistic}{the value of the standard deviate of Lee's L.}
  \item{p.value}{the p-value of the test.}
  \item{estimate}{the value of the observed Lee's L, its expectation and variance under the method assumption.}
  \item{alternative}{a character string describing the alternative hypothesis.}
  \item{method}{a character string giving the assumption used for calculating the standard deviate.}
  \item{data.name}{a character string giving the name(s) of the data.}
}
\note{See Lee (2004) for details on how the asymptotic expectation and
variance of Lee's L is computed. In particular, check Lee (2004), table 1, 
page 1690. 

This test may fail for large datasets as the computation of the 
asymptotic expectation and variance requires the use of dense matrices.}


\references{Lee (2004). A generalized significance testing method for global
measures of spatial association: an extension of the Mantel test.  Environment
and Planning A 2004, volume 36, pages 1687 - 1703}

\author{Roger Bivand and Virgilio Gómez-Rubio \email{Virgilio.Gomez@uclm.es}}

\seealso{\code{\link{lee}}, \code{\link{lee.mc}}, \code{\link{listw2U}}}

\examples{
data(oldcol)
col.W <- nb2listw(COL.nb, style="W")
crime <- COL.OLD$CRIME

lee.test(crime, crime, col.W, zero.policy=TRUE)

#Test with missing values
x<-crime
y<-crime
x[1:5]<-NA
y[3:7]<-NA

lee.test(x, y, col.W, zero.policy=TRUE, na.action=na.omit)
#  lee.test(x, y, col.W, zero.policy=TRUE)#Stops with an error



data(boston, package="spData")
lw<-nb2listw(boston.soi)

x<-boston.c$CMEDV
y<-boston.c$CRIM

lee.test(x, y, lw, zero.policy=TRUE, alternative="less")

#Test with missing values
x[1:5]<-NA
y[3:7]<-NA

lee.test(x, y, lw, zero.policy=TRUE, alternative="less", na.action=na.omit)














}
\keyword{spatial}
