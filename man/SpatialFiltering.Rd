\name{SpatialFiltering}
\alias{SpatialFiltering}
\alias{print.SFResult}
\alias{fitted.SFResult}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Semi-parametric spatial filtering}
\description{
  The function selects eigenvectors in a semi-parametric spatial filtering approach to removing spatial dependence from linear models. Selection is by brute force by finding the single eigenvector reducing the standard variate of Moran's I for regression residuals most, and continuing until no candidate eigenvector reduces the value by more than \code{tol}. It returns a summary table from the selection process and a matrix of selected eigenvectors for the specified model.
}
\usage{
SpatialFiltering(formula, lagformula, data = list(), nb, glist = NULL, style = "C",
 zero.policy = NULL, tol = 0.1, zerovalue = 1e-04, ExactEV = FALSE,
 symmetric = TRUE, alpha=NULL, alternative="two.sided", verbose=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a symbolic description of the model to be fit, assuming a spatial error representation; when lagformula is given, it should include only the response and the intercept term}
  \item{lagformula}{An extra one-sided formula to be used when a spatial lag representation is desired; the intercept is excluded within the function if present because it is part of the formula argument, but excluding it explicitly in the lagformula argument in the presence of factors generates a collinear model matrix}
  \item{data}{an optional data frame containing the variables in the model}
  \item{nb}{an object of class \code{nb}}
  \item{glist}{list of general weights corresponding to neighbours}
  \item{style}{\code{style} can take values W, B, C, U, and S}
  \item{zero.policy}{default NULL, use global option value; if FALSE stop with error for any empty neighbour sets, if TRUE permit the weights list to be formed with zero-length weights vectors}
  \item{tol}{tolerance value for convergence of spatial filtering}
  \item{zerovalue}{eigenvectors with eigenvalues of an absolute value smaller than zerovalue will be excluded in eigenvector search}
  \item{ExactEV}{Set ExactEV=TRUE to use exact expectations and variances rather than the expectation and variance of Moran's I from the previous iteration, default FALSE}
  \item{symmetric}{Should the spatial weights matrix be forced to symmetry, default TRUE}
  \item{alpha}{if not NULL, used instead of the tol= argument as a stopping rule to choose all eigenvectors up to and including the one with a probability value exceeding alpha.}
  \item{alternative}{a character string specifying the alternative hypothesis, must be one of greater, less or two.sided (default).}
  \item{verbose}{default NULL, use global option value; if TRUE report eigenvectors selected}
}
\value{
  An \code{SFResult} object, with:
  \item{selection}{a matrix summarising the selection of eigenvectors for inclusion, with columns:
    \describe{
      \item{Step}{Step counter of the selection procedure}
      \item{SelEvec}{number of selected eigenvector (sorted descending)}
      \item{Eval}{its associated eigenvalue}
      \item{MinMi}{value Moran's I for residual autocorrelation}
      \item{ZMinMi}{standardized value of Moran's I assuming a normal approximation}
      \item{pr(ZI)}{probability value of the permutation-based standardized deviate for the given value of the alternative argument}
      \item{R2}{R\^2 of the model including exogenous variables and eigenvectors}
      \item{gamma}{regression coefficient of selected eigenvector in fit}
    }
    The first row is the value at the start of the search
  }
  \item{dataset}{a matrix of the selected eigenvectors in order of selection}
}
\references{Tiefelsdorf M, Griffith DA. (2007) Semiparametric Filtering of Spatial Autocorrelation: The Eigenvector Approach. Environment and Planning A, 39 (5) 1193 - 1221.}
\author{Yongwan Chun, Michael Tiefelsdorf, Roger Bivand}

\seealso{\code{\link[stats]{lm}}, \code{\link[base]{eigen}}, \code{\link{nb2listw}}, \code{\link{listw2U}}}
\examples{
if (require(rgdal, quietly=TRUE)) {
example(columbus, package="spData")
lmbase <- lm(CRIME ~ INC + HOVAL, data=columbus)
sarcol <- SpatialFiltering(CRIME ~ INC + HOVAL, data=columbus,
 nb=col.gal.nb, style="W", ExactEV=TRUE)
sarcol
lmsar <- lm(CRIME ~ INC + HOVAL + fitted(sarcol), data=columbus)
lmsar
anova(lmbase, lmsar)
lm.morantest(lmsar, nb2listw(col.gal.nb))
lagcol <- SpatialFiltering(CRIME ~ 1, ~ INC + HOVAL - 1, data=columbus,
 nb=col.gal.nb, style="W")
lagcol
lmlag <- lm(CRIME ~ INC + HOVAL + fitted(lagcol), data=columbus)
lmlag
anova(lmbase, lmlag)
lm.morantest(lmlag, nb2listw(col.gal.nb))
}
}
\keyword{spatial}

