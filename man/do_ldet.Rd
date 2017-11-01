\name{do_ldet}
\alias{do_ldet}
\alias{jacobianSetup}
\alias{eigen_setup}
\alias{eigen_pre_setup}
\alias{mcdet_setup}
\alias{cheb_setup}
\alias{spam_setup}
\alias{spam_update_setup}
\alias{Matrix_setup}
\alias{Matrix_J_setup}
\alias{LU_setup}
\alias{LU_prepermutate_setup}
\alias{moments_setup}
\alias{SE_classic_setup}
\alias{SE_whichMin_setup}
\alias{SE_interp_setup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Spatial regression model Jacobian computations}

\description{These functions are made available in the package namespace for other developers, and are not intended for users. They provide a shared infrastructure for setting up data for Jacobian computation, and then for caclulating the Jacobian, either exactly or approximately, in maximum likelihood fitting of spatial regression models. The techniques used are the exact eigenvalue, Cholesky decompositions (Matrix, spam), and LU ones, with Chebyshev and Monte Carlo approximations; moments use the methods due to Martin and Smirnov/Anselin.}

\usage{
do_ldet(coef, env, which=1)
jacobianSetup(method, env, con, pre_eig=NULL, trs=NULL, interval=NULL, which=1)
cheb_setup(env, q=5, which=1)
mcdet_setup(env, p=16, m=30, which=1)
eigen_setup(env, which=1)
eigen_pre_setup(env, pre_eig, which=1)
spam_setup(env, pivot="MMD", which=1)
spam_update_setup(env, in_coef=0.1, pivot="MMD", which=1)
Matrix_setup(env, Imult, super=as.logical(NA), which=1)
Matrix_J_setup(env, super=FALSE, which=1)
LU_setup(env, which=1)
LU_prepermutate_setup(env, coef=0.1, order=FALSE, which=1)
moments_setup(env, trs=NULL, m, p, type="MC", correct=TRUE, trunc=TRUE, eq7=TRUE, which=1)
SE_classic_setup(env, SE_method="LU", p=16, m=30, nrho=200, interpn=2000,
 interval=c(-1,0.999), SElndet=NULL, which=1)
SE_whichMin_setup(env, SE_method="LU", p=16, m=30, nrho=200, interpn=2000,
 interval=c(-1,0.999), SElndet=NULL, which=1)
SE_interp_setup(env, SE_method="LU", p=16, m=30, nrho=200,
 interval=c(-1,0.999), which=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coef}{spatial coefficient value}
  \item{env}{environment containing pre-computed objects, fixed after assignment in setup functions}
  \item{which}{default 1; if 2, use second listw object}
  \item{method}{string value, used by \code{jacobianSetup} to choose method}
  \item{con}{control list passed from model fitting function and parsed in \code{jacobianSetup} to set environment variables for method-specific setup}
  \item{pre_eig}{pre-computed eigenvalues of length n}
  \item{q}{Chebyshev approximation order; default in calling spdep functions is 5, here it cannot be missing and does not have a default}
  \item{p}{Monte Carlo approximation number of random normal variables; default calling spdep functions is 16, here it cannot be missing and does not have a default}
  \item{m}{Monte Carlo approximation number of series terms; default in calling spdep functions is 30, here it cannot be missing and does not have a default; \code{m} serves the same purpose in the moments method}
  \item{pivot}{default \dQuote{MMD}, may also be \dQuote{RCM} for Cholesky decompisition using spam}
  \item{in_coef}{fill-in initiation coefficient value, default 0.1}
  \item{Imult}{see \code{\link[Matrix]{Cholesky}}; numeric scalar which defaults to zero. The matrix that is decomposed is A+m*I where m is the value of Imult and I is the identity matrix of order ncol(A). Default in calling spdep functions is 2, here it cannot be missing and does not have a default, but is rescaled for binary weights matrices in proportion to the maximim row sum in those calling functions}
  \item{super}{see \code{\link[Matrix]{Cholesky}}; logical scalar indicating is a supernodal decomposition should be created.  The alternative is a simplicial decomposition. Default in calling spdep functions is FALSE for \dQuote{Matrix_J} and \code{as.logical(NA)} for \dQuote{Matrix}.  Setting it to NA leaves the choice to a CHOLMOD-internal heuristic}
  \item{order}{default FALSE; used in LU_prepermutate, note warnings given for \code{lu} method}
  \item{trs}{A numeric vector of \code{m} traces, as from \code{trW}}
  \item{type}{moments trace type, see \code{\link{trW}}}
  \item{correct}{default TRUE: use Smirnov correction term, see \code{\link{trW}}}
  \item{trunc}{default TRUE: truncate Smirnov correction term, see \code{\link{trW}}}
  \item{eq7}{default TRUE}{use equation 7 in Smirnov and Anselin (2009), if FALSE no unit root correction}
  \item{SE_method}{default \dQuote{LU}, alternatively \dQuote{MC}; underlying lndet method to use for generating SE toolbox emulation grid}
  \item{nrho}{default 200, number of lndet values in first stage SE toolbox emulation grid}
  \item{interval}{default c(-1,0.999) if interval argument NULL, bounds for SE toolbox emulation grid}
  \item{interpn}{default 2000, number of lndet values to interpolate in second stage SE toolbox emulation grid}
  \item{SElndet}{default NULL, used to pass a pre-computed two-column matrix of coefficient values and corresponding interpolated lndet values}
}

\details{Since environments are containers in the R workspace passed by reference rather than by value, they are useful for passing objects to functions called in numerical optimisation, here for the maximum likelihood estimation of spatial regression models. This technique can save a little time on each function call, balanced against the need to access the objects in the environment inside the function. The environment should contain a \code{family} string object either \dQuote{SAR}, \dQuote{CAR} or \dQuote{SMA} (used in \code{do_ldet} to choose spatial moving average in \code{spautolm}, and these specific objects before calling the set-up functions:

\describe{
  \item{eigen}{Classical Ord eigenvalue computations - either:
  \describe{
    \item{listw}{A listw spatial weights object}
    \item{can.sim}{logical scalar: can the spatial weights be made symmetric by similarity}
    \item{verbose}{logical scalar: legacy report print control, for historical reasons only}
  } or:
  \describe{
    \item{pre_eig}{pre-computed eigenvalues}
  }
  and assigns to the environment:
  \describe{
    \item{eig}{a vector of eigenvalues}
    \item{eig.range}{the search interval for the spatial coefficient}
    \item{method}{string: \dQuote{eigen}}
  }
}
  \item{Matrix}{Sparse matrix pre-computed Cholesky decomposition with fast updating:
  \describe{
    \item{listw}{A listw spatial weights object}
    \item{can.sim}{logical scalar: can the spatial weights be made symmetric by similarity}
  } and assigns to the environment:
  \describe{
    \item{csrw}{sparse spatial weights matrix}
    \item{nW}{negative sparse spatial weights matrix}
    \item{pChol}{a \dQuote{CHMfactor} from factorising \code{csrw} with \code{\link[Matrix]{Cholesky}}}
    \item{nChol}{a \dQuote{CHMfactor} from factorising \code{nW} with \code{\link[Matrix]{Cholesky}}}
    \item{method}{string: \dQuote{Matrix}}
  }
}
  \item{Matrix_J}{Standard Cholesky decomposition without updating:
  \describe{
    \item{listw}{A listw spatial weights object}
    \item{can.sim}{logical scalar: can the spatial weights be made symmetric by similarity}
    \item{n}{number of spatial objects}
  } and assigns to the environment:
  \describe{
    \item{csrw}{sparse spatial weights matrix}
    \item{I}{sparse identity matrix}
    \item{super}{the value of the \code{super} argument}
    \item{method}{string: \dQuote{Matrix_J}}
  }
}
  \item{spam}{Standard Cholesky decomposition without updating:
  \describe{
    \item{listw}{A listw spatial weights object}
    \item{can.sim}{logical scalar: can the spatial weights be made symmetric by similarity}
    \item{n}{number of spatial objects}
  } and assigns to the environment:
  \describe{
    \item{csrw}{sparse spatial weights matrix}
    \item{I}{sparse identity matrix}
    \item{pivot}{string --- pivot method}
    \item{method}{string: \dQuote{spam}}
  }
}
  \item{spam_update}{Pre-computed Cholesky decomposition with updating:
  \describe{
    \item{listw}{A listw spatial weights object}
    \item{can.sim}{logical scalar: can the spatial weights be made symmetric by similarity}
    \item{n}{number of spatial objects}
  } and assigns to the environment:
  \describe{
    \item{csrw}{sparse spatial weights matrix}
    \item{I}{sparse identity matrix}
    \item{csrwchol}{A Cholesky decomposition for updating}
    \item{method}{string: \dQuote{spam}}
  }
}
  \item{LU}{Standard LU decomposition without updating:
  \describe{
    \item{listw}{A listw spatial weights object}
    \item{n}{number of spatial objects}
  } and assigns to the environment:
  \describe{
    \item{W}{sparse spatial weights matrix}
    \item{I}{sparse identity matrix}
    \item{method}{string: \dQuote{LU}}
  }
}
  \item{LU_prepermutate}{Standard LU decomposition with updating (pre-computed fill-reducing permutation):
  \describe{
    \item{listw}{A listw spatial weights object}
    \item{n}{number of spatial objects}
  } and assigns to the environment:
  \describe{
    \item{W}{sparse spatial weights matrix}
    \item{lu_order}{order argument to lu}
    \item{pq}{2-column matrix for row and column permutation for fill-reduction}
    \item{I}{sparse identity matrix}
    \item{method}{string: \dQuote{LU}}
  }
}
  \item{MC}{Monte Carlo approximation:
  \describe{
    \item{listw}{A listw spatial weights object}
  } and assigns to the environment:
  \describe{
    \item{clx}{list of Monte Carlo approximation terms  (the first two simulated traces are replaced by their analytical equivalents)}
    \item{W}{sparse spatial weights matrix}
    \item{method}{string: \dQuote{MC}}
  }
}
  \item{cheb}{Chebyshev approximation:
  \describe{
    \item{listw}{A listw spatial weights object}
  } and assigns to the environment:
  \describe{
    \item{trT}{vector of Chebyshev approximation terms}
    \item{W}{sparse spatial weights matrix}
    \item{method}{string: \dQuote{Chebyshev}}
  }
}
  \item{moments}{moments approximation:
  \describe{
    \item{listw}{A listw spatial weights object}
    \item{can.sim}{logical scalar: can the spatial weights be made symmetric by similarity}
  } and assigns to the environment:
  \describe{
    \item{trs}{vector of traces, possibly approximated}
    \item{q12}{integer vector of length 2, unit roots terms, ignored until 0.5-52}
    \item{eq7}{logical scalar: use equation 7}
    \item{correct}{logical scalar: use Smirnov correction term}
    \item{trunc}{logical scalar: truncate Smirnov correction term}
    \item{method}{string: \dQuote{moments}}
  }
}
  \item{SE_classic}{:
  \describe{
    \item{listw}{A listw spatial weights object}
    \item{n}{number of spatial objects}
  } and assigns to the environment:
  \describe{
    \item{detval}{two column matrix of lndet grid values}
    \item{method}{string: \dQuote{SE_classic}}
    \item{SE_method}{string: \dQuote{LU} or \dQuote{MC}}
  }
}
  \item{SE_whichMin}{:
  \describe{
    \item{listw}{A listw spatial weights object}
    \item{n}{number of spatial objects}
  } and assigns to the environment:
  \describe{
    \item{detval}{two column matrix of lndet grid values}
    \item{method}{string: \dQuote{SE_whichMin}}
    \item{SE_method}{string: \dQuote{LU} or \dQuote{MC}}
  }
}
  \item{SE_interp}{:
  \describe{
    \item{listw}{A listw spatial weights object}
    \item{n}{number of spatial objects}
  } and assigns to the environment:
  \describe{
    \item{fit}{fitted spline object from which to predict lndet values}
    \item{method}{string: \dQuote{SE_interp}}
    \item{SE_method}{string: \dQuote{LU} or \dQuote{MC}}
  }
}
}

Some set-up functions may also assign \code{similar} to the environment if the weights were made symmetric by similarity.

Three set-up functions emulate the behaviour of the Spatial Econometrics toolbox (March 2010) maximum likelihood lndet grid performance. The toolbox lndet functions compute a smaller number of lndet values for a grid of coefficient values (spacing 0.01), and then interpolate to a finer grid of values (spacing 0.001). \dQuote{SE_classic}, which is an implementation of the SE toolbox code, for example in f_sar.m, appears to have selected a row in the grid matrix one below the correct row when the candidate coefficient value was between 0.005 and 0.01-fuzz, always rounding the row index down. A possible alternative is to choose the index that is closest to the candidate coefficient value (\dQuote{SE_whichMin}). Another alternative is to fit a spline model to the first stage coarser grid, and pass this fitted model to the log likelihood function to make a point prediction using the candidate coefficient value, rather than finding the grid index (\dQuote{SE_interp}).

}

\value{\code{do_ldet} returns the value of the Jacobian for the calculation method recorded in the environment argument, and for the Monte Carlo approximation, returns a measure of the spread of the approximation as an \dQuote{sd} attribute; the remaining functions modify the environment in place as a side effect and return nothing.}

\references{LeSage J and RK Pace (2009) Introduction to Spatial Econometrics. CRC Press, Boca Raton, pp. 77--110.

Bivand, R. S., Hauke, J., and Kossowski, T. (2013). Computing the Jacobian in Gaussian spatial autoregressive models: An illustrated comparison of available methods. \emph{Geographical Analysis}, 45(2), 150-179.
}

\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{spautolm}}, \code{\link{lagsarlm}}, \code{\link{errorsarlm}}, \code{\link[Matrix]{Cholesky}}}

\examples{
data(boston, package="spData")
lw <- nb2listw(boston.soi)
can.sim <- spdep:::can.be.simmed(lw)
env <- new.env(parent=globalenv())
assign("listw", lw, envir=env)
assign("can.sim", can.sim, envir=env)
assign("similar", FALSE, envir=env)
assign("verbose", FALSE, envir=env)
assign("family", "SAR", envir=env)
eigen_setup(env)
get("similar", envir=env)
do_ldet(0.5, env)
rm(env)
env <- new.env(parent=globalenv())
assign("listw", lw, envir=env)
assign("can.sim", can.sim, envir=env)
assign("similar", FALSE, envir=env)
assign("verbose", FALSE, envir=env)
assign("family", "SAR", envir=env)
assign("n", length(boston.soi), envir=env)
eigen_pre_setup(env, pre_eig=eigenw(similar.listw(lw)))
do_ldet(0.5, env)
rm(env)
env <- new.env(parent=globalenv())
assign("listw", lw, envir=env)
assign("can.sim", can.sim, envir=env)
assign("similar", FALSE, envir=env)
assign("family", "SAR", envir=env)
assign("n", length(boston.soi), envir=env)
Matrix_setup(env, Imult=2, super=FALSE)
get("similar", envir=env)
do_ldet(0.5, env)
rm(env)
env <- new.env(parent=globalenv())
assign("listw", lw, envir=env)
assign("n", length(boston.soi), envir=env)
assign("can.sim", can.sim, envir=env)
assign("similar", FALSE, envir=env)
assign("family", "SAR", envir=env)
spam_setup(env)
get("similar", envir=env)
do_ldet(0.5, env)
rm(env)
env <- new.env(parent=globalenv())
assign("listw", lw, envir=env)
assign("n", length(boston.soi), envir=env)
assign("similar", FALSE, envir=env)
assign("family", "SAR", envir=env)
LU_setup(env)
get("similar", envir=env)
do_ldet(0.5, env)
rm(env)
env <- new.env(parent=globalenv())
assign("listw", lw, envir=env)
assign("n", length(boston.soi), envir=env)
assign("similar", FALSE, envir=env)
assign("family", "SAR", envir=env)
LU_prepermutate_setup(env)
get("similar", envir=env)
do_ldet(0.5, env)
rm(env)
env <- new.env(parent=globalenv())
assign("listw", lw, envir=env)
assign("similar", FALSE, envir=env)
assign("family", "SAR", envir=env)
cheb_setup(env, q=5)
get("similar", envir=env)
do_ldet(0.5, env)
rm(env)
env <- new.env(parent=globalenv())
assign("listw", lw, envir=env)
assign("n", length(boston.soi), envir=env)
assign("similar", FALSE, envir=env)
assign("family", "SAR", envir=env)
set.seed(12345)
mcdet_setup(env, p=16, m=30)
get("similar", envir=env)
do_ldet(0.5, env)
rm(env)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}

