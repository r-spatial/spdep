\name{lagmess}
\alias{lagmess}
\alias{print.lagmess}
\alias{print.summary.lagmess}
\alias{summary.lagmess}
\alias{residuals.lagmess}
\alias{deviance.lagmess}
\alias{coef.lagmess}
\alias{fitted.lagmess}
\alias{logLik.lagmess}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Matrix exponential spatial lag model}
\description{The function fits a matrix exponential spatial lag model, using \code{optim} to find the value of \code{alpha}, the spatial coefficient.}
\usage{
lagmess(formula, data = list(), listw, zero.policy = NULL, na.action = na.fail,
 q = 10, start = -2.5, control=list(), method="BFGS", verbose=NULL,
 use_expm=FALSE)
\method{summary}{lagmess}(object, ...)
\method{print}{lagmess}(x, ...)
\method{print}{summary.lagmess}(x, digits = max(5, .Options$digits - 3),
    signif.stars = FALSE, ...)
\method{residuals}{lagmess}(object, ...)
\method{deviance}{lagmess}(object, ...)
\method{coef}{lagmess}(object, ...)
\method{fitted}{lagmess}(object, ...)
\method{logLik}{lagmess}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a symbolic description of the model to be fit. The details 
of model specification are given for \code{lm()}}
  \item{data}{an optional data frame containing the variables in the model. 
By default the variables are taken from the environment which the function 
is called.}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{zero.policy}{default NULL, use global option value; if TRUE assign zero to the lagged value of zones without 
neighbours, if FALSE assign NA - causing \code{lagmess()} to terminate with an error}
  \item{na.action}{a function (default \code{options("na.action")}), can also be \code{na.omit} or \code{na.exclude} with consequences for residuals and fitted values - in these cases the weights list will be subsetted to remove NAs in the data. It may be necessary to set zero.policy to TRUE because this subsetting may create no-neighbour observations. Note that only weights lists created without using the glist argument to \code{nb2listw} may be subsetted.}
  \item{q}{default 10; number of powers of the spatial weights to use}
  \item{start}{starting value for numerical optimization, should be a small negative number}
  \item{control}{control parameters passed to \code{optim}}
  \item{method}{default \code{BFGS}, method passed to \code{optim}}
  \item{verbose}{default NULL, use global option value; if TRUE report function values during optimization}
  \item{use_expm}{default FALSE; if TRUE use \code{expm::expAtv} instead of a truncated power series of W}
  \item{x,object}{Objects of classes \code{lagmess} or \code{summary.lagmess} to be passed to methods}
  \item{digits}{the number of significant digits to use when printing}
  \item{signif.stars}{logical. If TRUE, "significance stars" are printed
          for each coefficient.}
  \item{\dots}{further arguments passed to or from other methods}
}

\details{The underlying spatial lag model:

\deqn{y = \rho W y + X \beta + \varepsilon}{y = rho W y + X beta + e}

where \eqn{\rho}{rho} is the spatial parameter may be fitted by maximum likelihood. In that case, the log likelihood function includes the logartithm of cumbersome Jacobian term \eqn{|I - \rho W|}{|I - rho W|}. If we rewrite the model as:

\deqn{S y = X \beta + \varepsilon}{S y = X beta + e}

we see that in the ML case \eqn{S y = (I - \rho W) y}{S y = (I - rho W) y}. If W is row-stochastic, S may be expressed as a linear combination of row-stochastic matrices. By pre-computing the matrix \eqn{[y Wy, W^2y, ..., W^{q-1}y]}{[y Wy, W^2y, ..., W^{q-1}y]}, the term \eqn{S y (\alpha)}{S y (alpha)} can readily be found by numerical optimization using the matrix exponential approach. \eqn{\alpha}{alpha} and \eqn{\rho}{rho} are related as \eqn{\rho = 1 - \exp{\alpha}}{rho = 1 - exp(alpha)}, conditional on the number of matrix power terms taken \code{q}.}

\value{
The function returns an object of class \code{lagmess} with components:
  \item{lmobj}{the \code{lm} object returned after fitting \code{alpha}}
  \item{alpha}{the spatial coefficient}
  \item{alphase}{the standard error of the spatial coefficient using the numerical Hessian}
  \item{rho}{the value of \code{rho} implied by \code{alpha}}
  \item{bestmess}{the object returned by \code{optim}}
  \item{q}{the number of powers of the spatial weights used}
  \item{start}{the starting value for numerical optimization used}
  \item{na.action}{(possibly) named vector of excluded or omitted observations if non-default na.action argument used}
  \item{nullLL}{the log likelihood of the aspatial model for the same data}
}

\references{J. P. LeSage and R. K. Pace (2007) A matrix exponential specification. Journal of Econometrics, 140, 190-214; J. P. LeSage and R. K. Pace (2009) Introduction to Spatial Econometrics. CRC Press, Chapter 9.}

\author{Roger Bivand \email{Roger.Bivand@nhh.no} and Eric Blankmeyer}

\seealso{\code{\link{lagsarlm}}, \code{\link[stats]{optim}}}

\examples{
data(baltimore, package="spData")
baltimore$AGE <- ifelse(baltimore$AGE < 1, 1, baltimore$AGE)
lw <- nb2listw(knn2nb(knearneigh(cbind(baltimore$X, baltimore$Y), k=7)))
obj1 <- lm(log(PRICE) ~ PATIO + log(AGE) + log(SQFT),
 data=baltimore)
lm.morantest(obj1, lw)
lm.LMtests(obj1, lw, test="all")
system.time(obj2 <- lagmess(log(PRICE) ~ PATIO + log(AGE) + log(SQFT), data=baltimore, listw=lw))
summary(obj2)
system.time(obj2a <- lagmess(log(PRICE) ~ PATIO + log(AGE) + log(SQFT), data=baltimore, listw=lw,
 use_expm=TRUE))
summary(obj2a)
obj3 <- lagsarlm(log(PRICE) ~ PATIO + log(AGE) + log(SQFT), data=baltimore, listw=lw)
summary(obj3)
data(boston, package="spData")
lw <- nb2listw(boston.soi)
gp2 <- lagsarlm(log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + I(RM^2)
 +  AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT),
 data=boston.c, lw, method="Matrix")
summary(gp2)
gp2a <- lagmess(CMEDV ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + I(RM^2)
 +  AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT),
 data=boston.c, lw)
summary(gp2a)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}

