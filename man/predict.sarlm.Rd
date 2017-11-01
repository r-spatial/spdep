% Copyright 2002 by Roger S. Bivand, 2015 Martin Gubri
\name{predict.sarlm}
\alias{predict.sarlm}
\alias{predict.SLX}
\alias{print.sarlm.pred}
\alias{as.data.frame.sarlm.pred}

\title{Prediction for spatial simultaneous autoregressive linear
model objects}
\description{
  \code{predict.sarlm()} calculates predictions as far as is at present possible for for spatial simultaneous autoregressive linear
model objects, using Haining's terminology for decomposition into
trend, signal, and noise, or other types of predictors --- see references.  
}
\usage{
\method{predict}{sarlm}(object, newdata = NULL, listw = NULL, pred.type = "TS", all.data = FALSE,
 zero.policy = NULL, legacy = TRUE, legacy.mixed = FALSE, power = NULL, order = 250,
 tol = .Machine$double.eps^(3/5), spChk = NULL, ...)
\method{predict}{SLX}(object, newdata, listw, zero.policy=NULL, ...)
\method{print}{sarlm.pred}(x, ...)
\method{as.data.frame}{sarlm.pred}(x, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{\code{sarlm} object returned by \code{lagsarlm}, 
\code{errorsarlm} or \code{sacsarlm}, the method for SLX objects takes the output of \code{lmSLX}}
  \item{newdata}{data frame in which to predict --- if NULL, predictions are
for the data on which the model was fitted. Should have row names corresponding to region.id. If row names are exactly the same than the ones used for training, it uses in-sample predictors for forecast.  See \sQuote{Details}}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}. In the out-of-sample prediction case (ie. if newdata is not NULL), if \code{legacy.mixed=FALSE} or if \code{pred.type!="TS"}, it should include both in-sample and out-of-sample spatial units. In this case, if regions of the listw are not in the correct order, they are reordered.  See \sQuote{Details}}
  \item{pred.type}{predictor type --- default \dQuote{TS}, use decomposition into  
trend, signal, and noise ; other types available depending on \code{newdata}. If \code{newdata=NULL} (in-sample prediction), \dQuote{TS}, \dQuote{trend}, \dQuote{TC} and \dQuote{BP} are available. If \code{newdata} is not NULL and its row names are the same than the \code{data} used to fit the model (forecast case), \dQuote{TS}, \dQuote{trend} and \dQuote{TC} are available. In other cases (out-of-sample prediction), \dQuote{TS}, \dQuote{trend}, \dQuote{KP1}, \dQuote{KP2}, \dQuote{KP3}, \dQuote{KP4}, \dQuote{KP5}, \dQuote{TC}, \dQuote{BP}, \dQuote{BPW}, \dQuote{BPN}, \dQuote{TS1}, \dQuote{TC1}, \dQuote{BP1}, \dQuote{BPW1} and \dQuote{BPN1} are available.  See \sQuote{Details} and references}
  \item{all.data}{(only applies to \code{pred.type="TC"} and newdata is not NULL) default FALSE: return predictions only for newdata units, if TRUE return predictions for all data units.  See \sQuote{Details}}
  \item{zero.policy}{default NULL, use global option value; if TRUE assign zero to the lagged value of zones without 
neighbours, if FALSE (default) assign NA - causing the function to 
terminate with an error}
  \item{legacy}{(only applies to lag and Durbin (mixed) models for \code{pred.type="TS"}) default TRUE: use ad-hoc predictor, if FALSE use DGP-based predictor}
  \item{legacy.mixed}{(only applies to mixed models if newdata is not NULL) default FALSE: compute lagged variables from both in-sample and out-of-sample units with \eqn{[W X]_O}{[WX]o} and \eqn{[W X]_S}{[WX]s} where \code{X=cbind(Xs, Xo)}, if TRUE compute lagged variables independantly between in-sample and out-of-sample units with \eqn{W_{OO} X_O}{Woo Xo} and \eqn{W_{SS} X_S}{Wss Xs} }
  \item{power}{(only applies to lag and Durbin (mixed) models for \dQuote{TS}, \dQuote{KP1}, \dQuote{KP2}, \dQuote{KP3}, \dQuote{TC}, \dQuote{TC1}, \dQuote{BP}, \dQuote{BP1}, \dQuote{BPN}, \dQuote{BPN1}, \dQuote{BPW} and \dQuote{BPW1} types) use \code{powerWeights}, if default NULL, set FALSE if \code{object$method} is \dQuote{eigen}, otherwise TRUE}
  \item{order}{power series maximum limit if \code{power} is TRUE}
  \item{tol}{tolerance for convergence of power series if \code{power} is TRUE}
  \item{spChk}{should the row names of data frames be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use \code{get.spChkOption()}}
  \item{x}{the object to be printed}
%   \item{pred.se}{logical, default FALSE, not available for error models}
%   \item{lagImpact}{object created by \code{impacts} with R sample draws from the fitted model}
  \item{...}{further arguments passed through}
}
\details{
The function supports three types of prediction.  In-sample prediction is the computation of predictors on the data used to fit the model (\code{newdata=NULL}).  Prevision, also called forecast, is the computation of some predictors (\dQuote{trend}, in-sample \dQuote{TC} and out-of-sample \dQuote{TS}) on the same spatial units than the ones used to fit the model, but with different observations of the variables in the model (row names of \code{newdata} should have the same row names than the data frame used to fit the model).  And out-of-sample prediction is the computation of predictors on other spatial units than the ones used to fit the model (\code{newdata} has different row names).  For extensive definitions, see Goulard et al. (2017). % add that warning if out-of-sample with newdata containing one in-sample spatial unit

\code{pred.type} of predictors are available according to the model of \code{object} an to the type of prediction.  In the two following tables, \dQuote{yes} means that the predictor can be used with the model, \dQuote{no} means that \code{predict.sarlm()} will stop with an error, and \dQuote{yes*} means that the predictor is not designed for the specified model, but it can be used with \code{predict.sarlm()}.  In the last case, be careful with the computation of a inappropriate predictor.

\emph{In-sample predictors by models}

\tabular{lccc}{
  pred.type \tab sem (mixed) \tab lag (mixed) \tab sac (mixed) \cr
        \tab     \tab  \tab  \cr
  \dQuote{trend} \tab  yes \tab  yes \tab  yes \cr
  \dQuote{TS}   \tab  yes \tab  yes \tab  no \cr
  \dQuote{TC}   \tab  no \tab  yes \tab  yes* \cr
  \dQuote{BP}   \tab  no \tab  yes \tab  yes* \cr
}

Note that only \dQuote{trend} and \dQuote{TC} are available for prevision.

\emph{Out-of-sample predictors by models}

\tabular{lccc}{
  pred.type \tab sem (mixed) \tab lag (mixed) \tab sac (mixed) \cr
        \tab     \tab  \tab  \cr
  \dQuote{trend} \tab  yes \tab  yes \tab  yes \cr
  \dQuote{TS}   \tab  yes \tab  yes \tab  no \cr
  \dQuote{TS1} or \dQuote{KP4}   \tab  no \tab  yes \tab  yes \cr
  \dQuote{TC}   \tab  no \tab  yes \tab  yes* \cr
  \dQuote{TC1} or \dQuote{KP1}  \tab  yes \tab  yes \tab  yes \cr
  \dQuote{BP}   \tab  no \tab  yes \tab  yes* \cr
  \dQuote{BP1}   \tab  no \tab  yes \tab  yes* \cr
  \dQuote{BPW}   \tab  no \tab  yes \tab  yes* \cr
  \dQuote{BPW1}   \tab  no \tab  yes \tab  yes* \cr
  \dQuote{BN}     \tab  no \tab  yes \tab  yes* \cr
  \dQuote{BPN1}   \tab  no \tab  yes \tab  yes* \cr
  \dQuote{KP2}   \tab  yes \tab  yes \tab  yes \cr
  \dQuote{KP3}   \tab  yes \tab  yes \tab  yes \cr
  \dQuote{KP5}   \tab  yes \tab  no \tab  yes* \cr
}

Values for \code{pred.type=} include \dQuote{TS1}, \dQuote{TC}, \dQuote{TC1}, \dQuote{BP}, \dQuote{BP1}, \dQuote{BPW}, \dQuote{BPW1}, \dQuote{BPN}, \dQuote{BPN1}, following the notation in Goulard et al. (2017), and for \code{pred.type=} \dQuote{KP1}, \dQuote{KP2}, \dQuote{KP3}, \dQuote{KP4}, \dQuote{KP5}, following the notation in Kelejian et al. (2007).  \code{pred.type="TS"} is described bellow and in Bivand (2002).  


In the following, the trend is the non-spatial smooth, the signal is the
spatial smooth, and the noise is the residual.  The fit returned by \code{pred.type="TS"} is the
sum of the trend and the signal.  

When \code{pred.type="TS"}, the function approaches prediction first by dividing invocations between 
those with or without newdata.  When no newdata is present, the response 
variable may be reconstructed as the sum of the trend, the signal, and the
noise (residuals).  Since the values of the response variable are known,
their spatial lags are used to calculate signal components (Cressie 1993, p. 564).  For the error
model, trend = \eqn{X \beta}{X beta}, and signal = \eqn{\lambda W y - 
\lambda W X \beta}{lambda W y - lambda W X beta}. For the lag and mixed
models, trend = \eqn{X \beta}{X beta}, and signal = \eqn{\rho W y}{rho W y}.

This approach differs from the design choices made in other software, for
example GeoDa, which does not use observations of the response variable,
and corresponds to the newdata situation described below.  

When however newdata is used for prediction, no observations of the response 
variable being predicted are available.  Consequently, while the trend
components are the same, the signal cannot take full account of the spatial
smooth.  In the error model and Durbin error model, the signal is set to zero, since the spatial smooth is expressed in terms of the error: 
\eqn{(I - \lambda W)^{-1} \varepsilon}{inv(I - lambda W) e}.

In the lag model, the signal can be expressed in the following way (for legacy=TRUE):

\deqn{(I - \rho W) y = X \beta + \varepsilon}{(I - rho W) y = X beta + e}
\deqn{y = (I - \rho W)^{-1} X \beta + (I - \rho W)^{-1} \varepsilon}{y = inv(I - rho W) X beta + inv(I - rho W) e}

giving a feasible signal component of:

\deqn{\rho W y = \rho W (I - \rho W)^{-1} X \beta}{rho W y = rho W inv(I - rho W) X beta}

For legacy=FALSE, the trend is computed first as:

\deqn{X \beta}{X beta}

next the prediction using the DGP:

\deqn{(I - \rho W)^{-1} X \beta}{inv(I - rho W) X beta}

and the signal is found as the difference between prediction and trend. The numerical results for the legacy and DGP methods are identical.  

%TODO: words missing?
setting the error term to zero.  This also means that predictions of the
signal component for lag and mixed models require the inversion of an 
n-by-n matrix.  

Because the outcomes of the spatial smooth on the error term are
unobservable, this means that the signal values for newdata are
incomplete.  In the mixed model, the spatially lagged RHS variables
influence both the trend and the signal, so that the root mean square
prediction error in the examples below for this case with newdata is
smallest, although the model was not the best fit.  

If \code{newdata} has more than one row, leave-one-out predictors (\code{pred.type=} include \dQuote{TS1}, \dQuote{TC1}, \dQuote{BP1}, \dQuote{BPW1}, \dQuote{BPN1}, \dQuote{KP1}, \dQuote{KP2}, \dQuote{KP3}, \dQuote{KP4}, \dQuote{KP5}) are computed separatly on each out-of-sample unit.  

\code{listw} should be provided except if \code{newdata=NULL} and \code{pred.type=} include \dQuote{TS}, \dQuote{trend}, or if \code{newdata} is not \code{NULL}, \code{pred.type="trend"} and \code{object} is not a mixed model.  % TODO: describe listw

\code{all.data} is useful when some out-of-sample predictors return different predictions for in-sample units, than the same predictor type computed only on in-sample data.  

}
\value{
  \code{predict.sarlm()} returns a vector of predictions with three attribute 
vectors of trend, signal (only for \code{pred.type="TS"}) and region.id values and two other attributes
of pred.type and call with class \code{sarlm.pred}. 

  \code{print.sarlm.pred()} is a print function for this class, printing and
returning a data frame with columns: "fit", "trend" and "signal" (when available) and with region.id as row names.
}

\references{Haining, R. 1990 \emph{Spatial data analysis in the social and environmental sciences}, Cambridge: Cambridge University Press, p. 258; Cressie, N. A. C. 1993 \emph{Statistics for spatial data}, Wiley, New York; Michel Goulard, Thibault Laurent & Christine Thomas-Agnan, 2017 \emph{About predictions in spatial autoregressive models: optimal and almost optimal strategies},  Spatial Economic Analysis Volume 12, Issue 2--3, 304--325, ; Kelejian, H. H. and Prucha, I. R. 2007 \emph{The relative efficiencies of various predictors in spatial econometric models containing spatial lags}, Regional Science and Urban Economics, Volume 37, Issue 3, 363--374; Bivand, R. 2002 \emph{Spatial econometrics functions in R: Classes and methods}, Journal of Geographical Systems, Volume 4, No. 4, 405--421}
\author{Roger Bivand \email{Roger.Bivand@nhh.no} and Martin Gubri}

\seealso{\code{\link{errorsarlm}}, \code{\link{lagsarlm}}, \code{\link{sacsarlm}}}

\examples{
data(oldcol)
lw <- nb2listw(COL.nb)
COL.lag.eig <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, lw)

COL.mix.eig <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, lw,
  type="mixed")
print(p1 <- predict(COL.mix.eig))
print(p2 <- predict(COL.mix.eig, newdata=COL.OLD, listw=lw, pred.type = "TS",
 legacy.mixed = TRUE))
AIC(COL.mix.eig)
sqrt(deviance(COL.mix.eig)/length(COL.nb))
sqrt(sum((COL.OLD$CRIME - as.vector(p1))^2)/length(COL.nb))
sqrt(sum((COL.OLD$CRIME - as.vector(p2))^2)/length(COL.nb))

COL.err.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, lw)
AIC(COL.err.eig)
sqrt(deviance(COL.err.eig)/length(COL.nb))
sqrt(sum((COL.OLD$CRIME - as.vector(predict(COL.err.eig)))^2)/length(COL.nb))
sqrt(sum((COL.OLD$CRIME - as.vector(predict(COL.err.eig, newdata=COL.OLD,
  listw=lw, pred.type = "TS")))^2)/length(COL.nb))

COL.SDerr.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, lw,
 etype="emixed")
AIC(COL.SDerr.eig)
sqrt(deviance(COL.SDerr.eig)/length(COL.nb))
sqrt(sum((COL.OLD$CRIME - as.vector(predict(COL.SDerr.eig)))^2)/length(COL.nb))
sqrt(sum((COL.OLD$CRIME - as.vector(predict(COL.SDerr.eig, newdata=COL.OLD,
  listw=lw, pred.type = "TS")))^2)/length(COL.nb))

AIC(COL.lag.eig)
sqrt(deviance(COL.lag.eig)/length(COL.nb))
sqrt(sum((COL.OLD$CRIME - as.vector(predict(COL.lag.eig)))^2)/length(COL.nb))
sqrt(sum((COL.OLD$CRIME - as.vector(predict(COL.lag.eig, newdata=COL.OLD,
  listw=lw, pred.type = "TS")))^2)/length(COL.nb))

p3 <- predict(COL.mix.eig, newdata=COL.OLD, listw=lw, pred.type = "TS",
 legacy=FALSE, legacy.mixed = TRUE)
all.equal(p2, p3, check.attributes=FALSE)
p4 <- predict(COL.mix.eig, newdata=COL.OLD, listw=lw, pred.type = "TS",
 legacy=FALSE, power=TRUE, legacy.mixed = TRUE)
all.equal(p2, p4, check.attributes=FALSE)
p5 <- predict(COL.mix.eig, newdata=COL.OLD, listw=lw, pred.type = "TS",
 legacy=TRUE, power=TRUE, legacy.mixed = TRUE)
all.equal(p2, p5, check.attributes=FALSE)

COL.SLX <- lmSLX(CRIME ~ INC + HOVAL, data=COL.OLD, listw=lw)
pslx0 <- predict(COL.SLX)
pslx1 <- predict(COL.SLX, newdata=COL.OLD, listw=lw)
all.equal(pslx0, pslx1)
COL.OLD1 <- COL.OLD
COL.OLD1$INC <- COL.OLD1$INC + 1
pslx2 <- predict(COL.SLX, newdata=COL.OLD1, listw=lw)
sum(coef(COL.SLX)[c(2,4)])
mean(pslx2-pslx1)
}


\keyword{spatial}
