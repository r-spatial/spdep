% Copyright 2004 by Roger S. Bivand
\name{p.adjustSP}
\alias{p.adjustSP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Adjust local association measures' p-values}
\description{
  Make an adjustment to local association measures' p-values based on the number of neighbours (+1) of each region, rather than the total number of regions.
}
\usage{
p.adjustSP(p, nb, method = "none")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{vector of p-values}
  \item{nb}{a list of neighbours of class \code{nb}}
  \item{method}{correction method as defined in \code{\link{p.adjust}}: "The adjustment methods include the Bonferroni correction ('"bonferroni"') in which the p-values are multiplied by the number of comparisons.  Four less conservative corrections are also included by Holm (1979) ('"holm"'), Hochberg (1988) ('"hochberg"'), Hommel (1988) ('"hommel"') and Benjamini & Hochberg (1995) ('"fdr"'), respectively. A pass-through option ('"none"') is also included."}
}

\value{
  A vector of corrected p-values using only the number of neighbours + 1.
}

\author{Danlin Yu and Roger Bivand \email{Roger.Bivand@nhh.no}}


\seealso{\code{\link{p.adjust}}, \code{\link{localG}}, \code{\link{localmoran}}}
\examples{
data(afcon, package="spData")
oid <- order(afcon$id)
resG <- as.vector(localG(afcon$totcon, nb2listw(include.self(paper.nb))))
non <- format.pval(pnorm(2*(abs(resG)), lower.tail=FALSE), 2)
bon <- format.pval(p.adjustSP(pnorm(2*(abs(resG)), lower.tail=FALSE),
 paper.nb, "bonferroni"), 2)
tot <- format.pval(p.adjust(pnorm(2*(abs(resG)), lower.tail=FALSE),
 "bonferroni", n=length(resG)), 2)
data.frame(resG, non, bon, tot, row.names=afcon$name)[oid,]
}
\keyword{spatial}
