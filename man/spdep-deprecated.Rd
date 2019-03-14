\name{spdep-deprecated}
\alias{spdep-deprecated}
\alias{lextrB}
\alias{lextrW}
\alias{lextrS}
\alias{l_max}
\alias{griffith_sone}
\alias{subgraph_eigenw}

\title{Deprecated Functions in Package \pkg{spdep}}

\description{
  These functions are provided for compatibility with older versions of
  \pkg{spdep} only, and may be defunct as soon as the next release. The
  functions have been moved to the \pkg{spreg} package.
}

\usage{
lextrB(lw, zero.policy = TRUE, control = list())
lextrW(lw, zero.policy=TRUE, control=list())
lextrS(lw, zero.policy=TRUE, control=list())
l_max(lw, zero.policy=TRUE, control=list())
griffith_sone(P, Q, type="rook")
subgraph_eigenw(nb, glist=NULL, style="W", zero.policy=NULL, quiet=NULL)
}

\arguments{
  \item{lw}{a binary symmetric \code{listw} object from, for example, \code{nb2listw} with style \dQuote{B} for \code{lextrB}, style \dQuote{W} for \code{lextrW} and style \dQuote{S} for \code{lextrS}; for \code{l_max}, the object may be asymmetric and does not have to be binary}
  \item{zero.policy}{default NULL, use global option value; if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA}
  \item{control}{a list of control arguments}
  \item{quiet}{default NULL, use global !verbose option value; set to FALSE for short summary}
  \item{P}{number of columns in the grid (number of units in a horizontal axis direction)}
  \item{Q}{number of rows in the grid (number of units in a vertical axis direction.)}
  \item{type}{\dQuote{rook} or \dQuote{queen}}
  \item{nb}{an object of class \code{nb}}
  \item{glist}{list of general weights corresponding to neighbours}
  \item{style}{\code{style} can take values \dQuote{W}, \dQuote{B}, \dQuote{C}, \dQuote{U}, \dQuote{minmax} and \dQuote{S}}

}
\section{Control arguments: lextrB/W/S, l_max}{
\describe{
  \item{trace}{report values in while loops, default NULL assuming FALSE; logical}
  \item{tol}{tolerance for breaking while loops, default \code{.Machine$double.eps^(1/2)}; numeric}
  \item{maxiter}{maximum number of iterations in while loops, default \code{6 * (length(lw$neighbours) - 2}; integer}
  \item{useC}{use C code, default TRUE, logical (not in \code{l_max})}
}
}

\details{
 Model-fitting functions and functions supporting model fitting are being moved to the \pkg{spreg} package.
The \code{griffith_sone} function function may be used, following Ord and Gasim (for references see Griffith and Sone (1995)), to calculate analytical eigenvalues for binary rook or queen contiguous neighbours where the data are arranged as a regular P times Q grid. The \code{subgraph_eigenw} function may be used when there are multiple graph components, of which the largest may be handled as a dense matrix. Here the eigenvalues are computed for each subgraph in turn, and catenated to reconstruct the complete set. The functions may be used to provide pre-computed eigenvalues for spatial regression functions.
}


\seealso{
  \code{\link[base]{Deprecated}}
}

\examples{
data(boston, package="spData")
ab.listb <- nb2listw(boston.soi, style="B")
er <- range(eigenw(ab.listb))
er
res_1 <- lextrB(ab.listb)
c(res_1)
run <- FALSE
if (require("RSpectra", quietly=TRUE)) run <- TRUE
if (run) {
B <- as(ab.listb, "CsparseMatrix")
eigs(B, k=1, which="SR")$values
}
if (run) {
eigs(B, k=1, which="LR")$values
}
k5 <- knn2nb(knearneigh(boston.utm, k=5))
c(l_max(nb2listw(k5, style="B")))
max(Re(eigenw(nb2listw(k5, style="B"))))
c(l_max(nb2listw(k5, style="C")))
max(Re(eigenw(nb2listw(k5, style="C"))))
ab.listw <- nb2listw(boston.soi, style="W")
er <- range(eigenw(similar.listw(ab.listw)))
er
res_1 <- lextrW(ab.listw)
c(res_1)
if (run) {
B <- as(similar.listw(ab.listw), "CsparseMatrix")
eigs(B, k=1, which="SR")$values
}
if (run) {
eigs(B, k=1, which="LR")$values
}
\dontrun{
ab.listw <- nb2listw(boston.soi, style="S")
er <- range(eigenw(similar.listw(ab.listw)))
er
res_1 <- lextrS(ab.listw)
c(res_1)
}
if (run) {
B <- as(similar.listw(ab.listw), "CsparseMatrix")
eigs(B, k=1, which="SR")$values
}
if (run) {
eigs(B, k=1, which="LR")$values
}
rg <- cell2nb(ncol=7, nrow=7, type="rook")
rg_eig <- eigenw(nb2listw(rg, style="B"))
rg_GS <- griffith_sone(P=7, Q=7, type="rook")
all.equal(rg_eig, rg_GS)
# subgraphs
data(oldcol)
crds <- cbind(COL.OLD$X, COL.OLD$Y)
k3 <- knn2nb(knearneigh(crds, k=3))
nc <- n.comp.nb(k3)
nc$nc
table(nc$comp.id)
k3eig <- eigenw(nb2listw(k3, style="W"))
k3eigSG <- subgraph_eigenw(k3, style="W")
all.equal(sort(k3eig), k3eigSG)

}

\keyword{spatial}

