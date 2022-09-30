\name{hotspot}
\alias{hotspot}
\alias{hotspot.default}
\alias{hotspot.localmoran}
\alias{hotspot.summary.localmoransad}
\alias{hotspot.data.frame.localmoranex}
\alias{hotspot.localG}
\alias{hotspot.localC}

\title{Cluster classifications for local indicators of spatial association}
\usage{
hotspot(obj, ...)

\method{hotspot}{default}(obj, ...)

\method{hotspot}{localmoran}(obj, Prname, cutoff=0.005, quadrant.type="mean",
 p.adjust="fdr", droplevels=TRUE, ...)
\method{hotspot}{summary.localmoransad}(obj, Prname, cutoff=0.005,
 quadrant.type="mean", p.adjust="fdr", droplevels=TRUE, ...)
\method{hotspot}{data.frame.localmoranex}(obj, Prname, cutoff=0.005,
 quadrant.type="mean", p.adjust="fdr", droplevels=TRUE, ...)

\method{hotspot}{localG}(obj, Prname, cutoff=0.005, p.adjust="fdr", droplevels=TRUE, ...)

\method{hotspot}{localC}(obj, Prname, cutoff=0.005, p.adjust="fdr", droplevels=TRUE, ...)
}
\arguments{
\item{obj}{An object of class \code{localmoran}, \code{localC} or \code{localG}}

\item{Prname}{A character string, the name of the column containing the probability values to be classified by cluster type if found \dQuote{interesting}}

\item{cutoff}{Default 0.005, the probability value cutoff larger than which the observation is not found \dQuote{interesting}}

\item{p.adjust}{Default \code{"fdr"}, the \code{p.adjust()} methood used, one of \code{c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")}}

\item{droplevels}{Default \code{TRUE}, should empty levels of the input cluster factor be dropped}

\item{quadrant.type}{Default \code{"mean"}, for \code{"localmoran"} objects only, can be \code{c("mean", "median", "pysal")} to partition the Moran scatterplot; \code{"mean"} partitions on the means of the variable and its spatial lag, \code{"median"} on medians of the variable and its spatial lag, \code{"pysal"} at zero for the centred variable and its spatial lag}

\item{...}{other arguments passed to methods.}

}
\description{
 Used to return a factor showing so-called cluster classification for local indicators of spatial association for local Moran's I, local Geary's C (and its multivariate variant) and local Getis-Ord G. This factor vector can be added to a spatial object for mapping.

}

\value{
 A factor showing so-called cluster classification for local indicators of spatial association.
}

\examples{
orig <- spData::africa.rook.nb
listw <- nb2listw(orig)
x <- spData::afcon$totcon

set.seed(1)
C <- localC_perm(x, listw)
Ch <- hotspot(C, Prname="Pr(z != E(Ci)) Sim", cutoff=0.05, p.adjust="none")
table(addNA(Ch))
set.seed(1)
I <- localmoran_perm(x, listw)
Ih <- hotspot(I, Prname="Pr(z != E(Ii)) Sim", cutoff=0.05, p.adjust="none")
table(addNA(Ih))
Is <- summary(localmoran.sad(lm(x ~ 1), nb=orig))
Ish <- hotspot(Is, Prname="Pr. (Sad)", cutoff=0.05, p.adjust="none")
table(addNA(Ish))
Ie <- as.data.frame(localmoran.exact(lm(x ~ 1), nb=orig))
Ieh <- hotspot(Ie, Prname="Pr. (exact)", cutoff=0.05, p.adjust="none")
table(addNA(Ieh))
set.seed(1)
G <- localG_perm(x, listw)
Gh <- hotspot(G, Prname="Pr(z != E(Gi)) Sim", cutoff=0.05, p.adjust="none")
table(addNA(Gh))
}

\author{
Roger Bivand
}
