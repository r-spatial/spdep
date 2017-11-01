\name{ME}
\alias{ME}
\alias{print.ME_res}
\alias{fitted.ME_res}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Moran eigenvector GLM filtering}
\description{
  The Moran eigenvector filtering function is intended to remove spatial autocorrelation from the residuals of generalised linear models. It uses brute force eigenvector selection to reach a subset of such vectors to be added to the RHS of the GLM model to reduce residual autocorrelation to below the specified alpha value. Since eigenvector selection only works on symmetric weights, the weights are made symmetric before the eigenvectors are found (from spdep 0.5-50).
}
\usage{
ME(formula, data, family = gaussian, weights, offset, listw,
 alpha=0.05, nsim=99, verbose=NULL, stdev=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a symbolic description of the model to be fit}
  \item{data}{an optional data frame containing the variables in the model}
  \item{family}{a description of the error distribution and link function to
          be used in the model}
  \item{weights}{an optional vector of weights to be used in the fitting process}
  \item{offset}{this can be used to specify an a priori known component to be included in the linear predictor during fitting}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{alpha}{used as a stopping rule to choose all eigenvectors up to and including the one with a p-value exceeding alpha}
  \item{nsim}{number of permutations for permutation bootstrap for finding p-values}
  \item{verbose}{default NULL, use global option value; if TRUE report eigenvectors selected}
  \item{stdev}{if TRUE, p-value calculated from bootstrap permutation standard deviate using \code{pnorm} with alternative="greater", if FALSE the Hope-type p-value}
}
\details{
  The eigenvectors for inclusion are chosen by calculating the empirical Moran's I values for the initial model plus each of the doubly centred symmetric spatial weights matrix eigenvectors in turn. Then the first eigenvector is chosen as that with the lowest Moran's I value. The procedure is repeated until the lowest remaining Moran's I value has a permutation-based probability value above alpha. The probability value is either Hope-type or based on using the mean and standard deviation of the permutations to calculate ZI based on the stdev argument.
}
\value{
An object of class \code{ME_res}:
  \item{selection}{a matrix summarising the selection of eigenvectors for inclusion, with columns:
    \describe{
      \item{Eigenvector}{number of selected eigenvector}
      \item{ZI}{permutation-based standardized deviate of Moran's I if stdev=TRUE}
      \item{pr(ZI)}{probability value: if stdev=TRUE of the permutation-based standardized deviate, if FALSE the Hope-type probability value, in both cases on-sided}
    }
    The first row is the value at the start of the search
  }
  \item{vectors}{a matrix of the selected eigenvectors in order of selection}
}
\references{Dray S, Legendre P and Peres-Neto PR (2005) Spatial modeling: a comprehensive framework for principle coordinate analysis of neigbbor matrices (PCNM), Ecological Modelling; Griffith DA and Peres-Neto PR (2006) Spatial modeling in ecology: the flexibility of eigenfunction spatial analyses.}
\author{Roger Bivand and Pedro Peres-Neto}

\seealso{\code{\link{SpatialFiltering}}, \code{\link[stats]{glm}}}

\examples{
\dontrun{
if (require(rgdal, quietly=TRUE)) {
example(columbus, package="spData")
lmbase <- lm(CRIME ~ INC + HOVAL, data=columbus)
lagcol <- SpatialFiltering(CRIME ~ 1, ~ INC + HOVAL, data=columbus,
 nb=col.gal.nb, style="W", alpha=0.1, verbose=TRUE)
lagcol
lmlag <- lm(CRIME ~ INC + HOVAL + fitted(lagcol), data=columbus)
anova(lmlag)
anova(lmbase, lmlag)
set.seed(123)
lagcol1 <- ME(CRIME ~ INC + HOVAL, data=columbus, family="gaussian",
 listw=nb2listw(col.gal.nb), alpha=0.1, verbose=TRUE)
lagcol1
lmlag1 <- lm(CRIME ~ INC + HOVAL + fitted(lagcol1), data=columbus)
anova(lmlag1)
anova(lmbase, lmlag1)
set.seed(123)
lagcol2 <- ME(CRIME ~ INC + HOVAL, data=columbus, family="gaussian",
 listw=nb2listw(col.gal.nb), alpha=0.1, stdev=TRUE, verbose=TRUE)
lagcol2
lmlag2 <- lm(CRIME ~ INC + HOVAL + fitted(lagcol2), data=columbus)
anova(lmlag2)
anova(lmbase, lmlag2)
example(nc.sids, package="spData")
glmbase <- glm(SID74 ~ 1, data=nc.sids, offset=log(BIR74),
 family="poisson")
set.seed(123)
MEpois1 <- ME(SID74 ~ 1, data=nc.sids, offset=log(BIR74),
 family="poisson", listw=nb2listw(ncCR85_nb, style="B"), alpha=0.2, verbose=TRUE)
MEpois1
glmME <- glm(SID74 ~ 1 + fitted(MEpois1), data=nc.sids, offset=log(BIR74),
 family="poisson")
anova(glmME, test="Chisq")
anova(glmbase, glmME, test="Chisq")
}
data(hopkins, package="spData")
hopkins_part <- hopkins[21:36,36:21]
hopkins_part[which(hopkins_part > 0, arr.ind=TRUE)] <- 1
hopkins.rook.nb <- cell2nb(16, 16, type="rook")
glmbase <- glm(c(hopkins_part) ~ 1, family="binomial")
set.seed(123)
MEbinom1 <- ME(c(hopkins_part) ~ 1, family="binomial",
 listw=nb2listw(hopkins.rook.nb, style="B"), alpha=0.2, verbose=TRUE)
glmME <- glm(c(hopkins_part) ~ 1 + fitted(MEbinom1), family="binomial")
anova(glmME, test="Chisq")
anova(glmbase, glmME, test="Chisq")
}
}
\keyword{spatial}

