\name{normalized.moran}
\alias{normalized.moran}
\title{Compute Normalized Moran's \eqn{I}}
\usage{
normalized.moran(
  x,
  listw,
  zero.policy = attr(listw, "zero.policy"),
  NAOK = FALSE
)
}
\arguments{
\item{x}{a numeric vector the same length as the neighbours list in listw}

\item{listw}{a \code{listw} object of style 'W' created for example by \code{nb2listw}}

\item{zero.policy}{default \code{attr(listw, "zero.policy")} as set when \code{listw} was created, if attribute not set, use global option value; if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA}

\item{NAOK}{if 'TRUE' then any 'NA' or 'NaN' or 'Inf' values in x are passed on to the foreign function. If 'FALSE', the presence of 'NA' or 'NaN' or 'Inf' values is regarded as an error.}
}
\value{
The normalized "\eqn{I_6}" Moran's \eqn{I} of the values in \code{x} is computed, as a numeric vector of length one.
}
\description{
A simple function to compute the normalized "\eqn{I_6}" Moran's \eqn{I} statistic proposed by Tillé et al. (2026), called by \code{normalized.moran.test} and \code{normalized.moran.mc}.
}
\section{Author}{
 Diego Giuliani \email{diego.giuliani@unitn.it}
}

\section{References}{
 Tillé Y., Giuliani D., Dickson M.M., Espa G. (2026) Revisiting the Normalization of the Moran's \eqn{I} index: A Correlation-based Approach with Inference, \emph{International Statistical Review}
}

\examples{
library(spdep)
data(oldcol)
col.W <- nb2listw(COL.nb, style="W")
crime <- COL.OLD$CRIME
normalized.moran(crime, col.W)
is.na(crime) <- sample(1:length(crime), 10)
normalized.moran(crime, col.W, NAOK=TRUE)
}
\seealso{
\link{normalized.moran.test}, \link{normalized.moran.mc}
}
