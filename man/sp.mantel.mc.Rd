\name{sp.mantel.mc}
\alias{sp.mantel.mc}
\alias{plot.mc.sim}

\title{Mantel-Hubert spatial general cross product statistic}
\description{
  A permutation test for the spatial general cross product statistic with
Moran (\eqn{C_{ij} = z_i z_j}{C(i,j) = z(i)*z(j)}), Geary 
(\eqn{C_{ij} = (z_i - z_j)^2}{C(i,j) = (z(i) - z(j))^2}), and Sokal 
(\eqn{C_{ij} = |z_i - z_j|}{C(i,j) = |z(i) - z(j)|}) criteria, for
\eqn{z_i = (x_i - \bar{x}) / \sigma_{x}}{z(i) = (x(i) - mean(x))/sd(x)}.
\code{plot.mc.sim} is a helper function to plot the outcomes of the
permutation test.
}
\usage{
sp.mantel.mc(var, listw, nsim, type = "moran", zero.policy = NULL,
 alternative = "greater", spChk=NULL, return_boot=FALSE)
\method{plot}{mc.sim}(x, xlim, xlab, main, sub, ..., ptype="density")
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{var}{a numeric vector the same length as the neighbours list in listw}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{nsim}{number of permutations}
  \item{type}{"moran", "geary" or "sokal" criteria for similarity}
  \item{zero.policy}{default NULL, use global option value; if TRUE assign zero to the lagged value of zones without 
  neighbours, if FALSE assign NA}
  \item{alternative}{a character string specifying the alternative hypothesis, 
  must be one of "greater" (default), or "less".}
  \item{spChk}{should the data vector names be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use \code{get.spChkOption()}}
  \item{return_boot}{return an object of class \code{boot} from the equivalent permutation bootstrap rather than an object of class \code{htest}}
  \item{x}{the object to be plotted}
  \item{xlim}{the range of the x axis}
  \item{xlab}{a title for the x axis}
  \item{main}{an overall title for the plot}
  \item{sub}{a sub title for the plot}
  \item{ptype}{either "density" or "hist"}
  \item{...}{further arguments passed through}
}

\value{
A list with class \code{htest} and \code{mc.sim} containing the following components:
  \item{statistic}{the value of the observed Geary's C.}
  \item{parameter}{the rank of the observed Geary's C.}
  \item{alternative}{a character string describing the alternative hypothesis.}
  \item{method}{a character string giving the method used.}
  \item{data.name}{a character string giving the name(s) of the data, and the 
  number of simulations.}
  \item{p.value}{the pseudo p-value of the test.}
  \item{res}{nsim simulated values of statistic, final value is observed 
  statistic}
  \item{estimate}{the mean and variance of the simulated distribution.}
}
\references{Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 22-24,
Haining, R. 1990 \emph{Spatial data analysis in the social
and environmental sciences}, Cambridge: Cambridge University Press, p. 230--1.
The function has been checked against general matrix code posted to the 
r-help list by Ben Bolker on 1 May 2001; another \code{mantel()} function
is in the vegan package.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{moran.mc}}, \code{\link{joincount.mc}}, 
\code{\link{geary.mc}}}

\examples{
data(oldcol)
sim1 <- sp.mantel.mc(COL.OLD$CRIME, nb2listw(COL.nb),
 nsim=99, type="geary", alternative="less")
sim1
plot(sim1)
sp.mantel.mc(COL.OLD$CRIME, nb2listw(COL.nb), nsim=99,
 type="sokal", alternative="less")
sp.mantel.mc(COL.OLD$CRIME, nb2listw(COL.nb), nsim=99,
 type="moran")
}
\keyword{spatial}
