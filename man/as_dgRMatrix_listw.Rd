\name{as_dgRMatrix_listw}
\alias{as_dgRMatrix_listw}
\alias{as_dsTMatrix_listw}
\alias{as_dsCMatrix_I}
\alias{as_dsCMatrix_IrW}
\alias{Jacobian_W}
\alias{coerce,listw,CsparseMatrix-method}
\alias{coerce,listw,RsparseMatrix-method}
\alias{coerce,listw,symmetricMatrix-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Interface between Matrix class objects and weights lists}
\description{
  Interface between Matrix class objects and weights lists
}
\usage{
as_dgRMatrix_listw(listw)
as_dsTMatrix_listw(listw)
as_dsCMatrix_I(n)
as_dsCMatrix_IrW(W, rho)
Jacobian_W(W, rho)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{W}{a \code{dsTMatrix} object created using \code{as_dsTMatrix_listw} from a symmetric \code{listw} object}
  \item{rho}{spatial regression coefficient}
  \item{n}{length of diagonal for identity matrix}
}

\value{
  Matrix package class objects
}

\author{Roger Bivand}

\examples{
example(NY_data)
W_C <- as(listw_NY, "CsparseMatrix")
W_R <- as(listw_NY, "RsparseMatrix")
W_S <- as(listw_NY, "symmetricMatrix")
n <- nrow(W_S)
I <- Diagonal(n)
rho <- 0.1
c(determinant(I - rho * W_S, logarithm=TRUE)$modulus)
sum(log(1 - rho * eigenw(listw_NY)))
nW <- - W_S
nChol <- Cholesky(nW, Imult=8)
n * log(rho) + (2 * c(determinant(update(nChol, nW, 1/rho))$modulus))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}

