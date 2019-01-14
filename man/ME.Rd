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
ME(formula, data, family = gaussian, weights, offset, na.action, listw,
 alpha=0.05, nsim=99, verbose=NULL, stdev=FALSE, zero.policy = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a symbolic description of the model to be fit}
  \item{data}{an optional data frame containing the variables in the model}
  \item{family}{a description of the error distribution and link function to
          be used in the model}
  \item{weights}{an optional vector of weights to be used in the fitting process}
  \item{offset}{this can be used to specify an a priori known component to be included in the linear predictor during fitting}
  \item{na.action}{a function (default \code{options("na.action")}), can also be \code{na.omit} or \code{na.exclude} with consequences for residuals and fitted values - in these cases the spatial weights list will be subsetted to remove NAs in the data. It may be necessary to set zero.policy to TRUE because this subsetting may create no-neighbour observations. Note that only weights lists created without using the glist argument to \code{nb2listw} may be subsetted.}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{alpha}{used as a stopping rule to choose all eigenvectors up to and including the one with a p-value exceeding alpha}
  \item{nsim}{number of permutations for permutation bootstrap for finding p-values}
  \item{verbose}{default NULL, use global option value; if TRUE report eigenvectors selected}
  \item{stdev}{if TRUE, p-value calculated from bootstrap permutation standard deviate using \code{pnorm} with alternative="greater", if FALSE the Hope-type p-value}
  \item{zero.policy}{default NULL, use global option value; if FALSE stop with error for any empty neighbour sets, if TRUE permit the weights list to be formed with zero-length weights vectors}
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
run <- FALSE
if (require(sf, quietly=TRUE)) run <- TRUE
if (run) {
columbus <- st_read(system.file("shapes/columbus.shp", package="spData")[1], quiet=TRUE)
col.gal.nb <- read.gal(system.file("weights/columbus.gal", package="spData")[1])
lmbase <- lm(CRIME ~ INC + HOVAL, data=columbus)
}
if (run) {
lagcol <- SpatialFiltering(CRIME ~ 1, ~ INC + HOVAL, data=columbus,
 nb=col.gal.nb, style="W", alpha=0.1, verbose=TRUE)
lagcol
}
if (run) {
lmlag <- lm(CRIME ~ INC + HOVAL + fitted(lagcol), data=columbus)
anova(lmbase, lmlag)
}
if (run) {
set.seed(123)
lagcol1 <- ME(CRIME ~ INC + HOVAL, data=columbus, family="gaussian",
 listw=nb2listw(col.gal.nb), alpha=0.1, verbose=TRUE)
lagcol1
}
if (run) {
lmlag1 <- lm(CRIME ~ INC + HOVAL + fitted(lagcol1), data=columbus)
anova(lmbase, lmlag1)
}
if (run) {
set.seed(123)
lagcol2 <- ME(CRIME ~ INC + HOVAL, data=columbus, family="gaussian",
 listw=nb2listw(col.gal.nb), alpha=0.1, stdev=TRUE, verbose=TRUE)
lagcol2
}
if (run) {
lmlag2 <- lm(CRIME ~ INC + HOVAL + fitted(lagcol2), data=columbus)
anova(lmbase, lmlag2)
}
if (run) {
NA.columbus <- columbus
NA.columbus$CRIME[20:25] <- NA
COL.ME.NA <- ME(CRIME ~ INC + HOVAL, data=NA.columbus, family="gaussian",
 listw=nb2listw(col.gal.nb), alpha=0.1, stdev=TRUE, verbose=TRUE,
 na.action=na.exclude)
COL.ME.NA$na.action
}
if (run) {
summary(lm(CRIME ~ INC + HOVAL + fitted(COL.ME.NA), data=NA.columbus,
 na.action=na.exclude))
}
if (run) {
nc.sids <- st_read(system.file("shapes/sids.shp", package="spData")[1], quiet=TRUE)
rn <- as.character(nc.sids$FIPS)
ncCC89_nb <- read.gal(system.file("weights/ncCC89.gal", package="spData")[1],
 region.id=rn)
ncCR85_nb <- read.gal(system.file("weights/ncCR85.gal", package="spData")[1],
 region.id=rn)
glmbase <- glm(SID74 ~ 1, data=nc.sids, offset=log(BIR74),
 family="poisson")
}
if (run) {
set.seed(123)
MEpois1 <- ME(SID74 ~ 1, data=nc.sids, offset=log(BIR74),
 family="poisson", listw=nb2listw(ncCR85_nb, style="B"), alpha=0.2, verbose=TRUE)
MEpois1
}
if (run) {
glmME <- glm(SID74 ~ 1 + fitted(MEpois1), data=nc.sids, offset=log(BIR74),
 family="poisson")
anova(glmME, test="Chisq")
}
if (run) {
anova(glmbase, glmME, test="Chisq")
}
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
\keyword{spatial}

