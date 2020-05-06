\name{moran.plot}
\alias{moran.plot}
\title{Moran scatterplot}
\description{
 A plot of spatial data against its spatially lagged values, augmented by
reporting the summary of influence measures for the linear relationship
between the data and the lag. If zero policy is TRUE, such observations
are also marked if they occur.
}
\usage{
moran.plot(x, listw, zero.policy=NULL, spChk=NULL, labels=NULL,
 xlab=NULL, ylab=NULL, quiet=NULL, plot=TRUE, return_df=TRUE, ...)
}
\arguments{
  \item{x}{a numeric vector the same length as the neighbours list in listw}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{zero.policy}{default NULL, use global option value; if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA}
  \item{spChk}{should the data vector names be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use \code{get.spChkOption()}}
  \item{labels}{character labels for points with high influence measures, if set to FALSE, no labels are plotted for points with large influence}
  \item{xlab}{label for x axis}
  \item{ylab}{label for x axis}
  \item{quiet}{default NULL, use !verbose global option value; if TRUE, output of summary of influence object suppressed}
  \item{plot}{default TRUE, if false, plotting is suppressed}
  \item{return_df}{default TRUE, invisibly return a data.frame object; if false invisibly return an influence measures object}
  \item{\dots}{further graphical parameters as in \code{par(..)}}
}

\value{
 The function returns a data.frame object with coordinates and influence measures if \code{return_df} is TRUE, or an influence object from \code{influence.measures}.
}
\references{Anselin, L. 1996. The Moran scatterplot as an ESDA tool to
assess local instability in spatial association. pp. 111--125 in
M. M. Fischer, H. J. Scholten and D. Unwin (eds) Spatial analytical
perspectives on GIS, London, Taylor and Francis; Anselin, L. 1995.
Local indicators of spatial association,
Geographical Analysis, 27, 93--115}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{localmoran}}, \code{\link[stats]{influence.measures}}}

\examples{
data(afcon, package="spData")
mp <- moran.plot(afcon$totcon, nb2listw(paper.nb),
 labels=as.character(afcon$name), pch=19)
moran.plot(as.vector(scale(afcon$totcon)), nb2listw(paper.nb),
 labels=as.character(afcon$name), xlim=c(-2, 4), ylim=c(-2,4), pch=19)
if (require(ggplot2, quietly=TRUE)) {
  xname <- attr(mp, "xname")
  ggplot(mp, aes(x=x, y=wx)) + geom_point(shape=1) + 
    geom_smooth(formula=y ~ x, method="lm") + 
    geom_hline(yintercept=mean(mp$wx), lty=2) + 
    geom_vline(xintercept=mean(mp$x), lty=2) + theme_minimal() + 
    geom_point(data=mp[mp$is_inf,], aes(x=x, y=wx), shape=9) +
    geom_text(data=mp[mp$is_inf,], aes(x=x, y=wx, label=labels, vjust=1.5)) +
    xlab(xname) + ylab(paste0("Spatially lagged ", xname))
}
}
\keyword{spatial}
