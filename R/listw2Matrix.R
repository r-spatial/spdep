# Copyright 2006-14 by Roger Bivand
#

setAs("listw", "CsparseMatrix", function(from) {
if (requireNamespace("spatialreg", quietly=TRUE)) {
  as(spatialreg::as_dgRMatrix_listw(from), "CsparseMatrix")
} else {
  as(as_dgRMatrix_listw(from), "CsparseMatrix")
#  stop("install the spatialreg package")
}
})
setAs("listw", "RsparseMatrix", function(from) {
if (requireNamespace("spatialreg", quietly=TRUE)) {
  spatialreg::as_dgRMatrix_listw(from)
} else {
  as_dgRMatrix_listw(from)
#  stop("install the spatialreg package")
}
})
setAs("listw", "symmetricMatrix", function(from) {
if (requireNamespace("spatialreg", quietly=TRUE)) {
  spatialreg::as_dsTMatrix_listw(from)
} else {
  as_dsTMatrix_listw(from)
#  stop("install the spatialreg package")
}
})


as_dgRMatrix_listw <- function(listw) {
    .Deprecated("spatialreg::as_dgRMatrix_listw", msg="Function as_dgRMatrix_listw moved to the spatialreg package")
#    if (!requireNamespace("spatialreg", quietly=TRUE))
#      stop("install the spatialreg package")
    if (requireNamespace("spatialreg", quietly=TRUE)) {
      return(spatialreg::as_dgRMatrix_listw(listw=listw))
    }
    warning("install the spatialreg package")
#  if (FALSE) {
	if(!inherits(listw, "listw")) stop("not a listw object")
	n <- length(listw$neighbours)
	cardw <- card(listw$neighbours)
	p0 <- as.integer(c(0, cumsum(cardw)))
	scard <- sum(cardw)
	z <- .Call("listw2dgR", listw$neighbours, listw$weights,
		as.integer(cardw), as.integer(scard), PACKAGE="spdep")
	res <- new("dgRMatrix", j=z[[1]], p=p0, Dim=as.integer(c(n, n)),
		x=z[[2]])
        colnames(res) <- attr(listw$neighbours, "region.id")
        rownames(res) <- colnames(res)
	res
}
#}

as_dsTMatrix_listw <- function(listw) {
    .Deprecated("spatialreg::as_dsTMatrix_listw", msg="Function as_dsTMatrix_listw moved to the spatialreg package")
#    if (!requireNamespace("spatialreg", quietly=TRUE))
#      stop("install the spatialreg package")
    if (requireNamespace("spatialreg", quietly=TRUE)) {
      return(spatialreg::as_dsTMatrix_listw(listw=listw))
    }
    warning("install the spatialreg package")
#  if (FALSE) {
	if (!inherits(listw, "listw")) stop("not a listw object")
	if (!is.symmetric.glist(listw$neighbours, listw$weights))
		stop("not a symmetric matrix")
	n <- length(listw$neighbours)
	cardw <- card(listw$neighbours)
	scard <- sum(cardw)
	if (scard %% 2 != 0) stop("odd neighbours sum")
	z <- .Call("listw2dsT", listw$neighbours, listw$weights,
		as.integer(cardw), as.integer(scard/2), PACKAGE="spdep")

	res <- new("dsTMatrix", i=z[[1]], j=z[[2]], Dim=as.integer(c(n, n)),
		x=z[[3]])
        colnames(res) <- attr(listw$neighbours, "region.id")
        rownames(res) <- colnames(res)
	res
}
#}

as_dsCMatrix_I <- function(n) {
    .Deprecated("spatialreg::as_dsCMatrix_I", msg="Function as_dsCMatrix_I moved to the spatialreg package")
#    if (!requireNamespace("spatialreg", quietly=TRUE))
#      stop("install the spatialreg package")
    if (requireNamespace("spatialreg", quietly=TRUE)) {
      return(spatialreg::as_dsCMatrix_I(n=n))
    }
    warning("install the spatialreg package")
#  if (FALSE) {
	if (n < 1) stop("matrix must have positive dimensions")
	as(as(Diagonal(n), "symmetricMatrix"), "CsparseMatrix")
}
#}

as_dsCMatrix_IrW <- function(W, rho) {
    .Deprecated("spatialreg::as_dsCMatrix_IrW", msg="Function as_dsCMatrix_IrW moved to the spatialreg package")
#    if (!requireNamespace("spatialreg", quietly=TRUE))
#      stop("install the spatialreg package")
    if (requireNamespace("spatialreg", quietly=TRUE)) {
      return(spatialreg::as_dsCMatrix_IrW(W=W, rho=rho))
    }
    warning("install the spatialreg package")
#  if (FALSE) {
	stopifnot(is(W, "symmetricMatrix"))
	n <- dim(W)[1]
	as_dsCMatrix_I(n) - rho * W
}
#}

Jacobian_W <- function(W, rho) {
    .Deprecated("spatialreg::Jacobian_W", msg="Function Jacobian_W moved to the spatialreg package")
#    if (!requireNamespace("spatialreg", quietly=TRUE))
#      stop("install the spatialreg package")
    if (requireNamespace("spatialreg", quietly=TRUE)) {
      return(spatialreg::Jacobian_W(W=W, rho=rho))
    }
    warning("install the spatialreg package")
#  if (FALSE) {
	sum(2*log(diag(chol(as_dsCMatrix_IrW(W, rho)))))
}
#}


listw2U_Matrix <- function(lw) {
    .Deprecated("spatialreg::listw2U_Matrix", msg="Function listw2U_Matrix moved to the spatialreg package")
#    if (!requireNamespace("spatialreg", quietly=TRUE))
#      stop("install the spatialreg package")
    if (requireNamespace("spatialreg", quietly=TRUE)) {
      return(spatialreg::listw2U_Matrix(lw=lw))
    }
    warning("install the spatialreg package")
#  if (FALSE) {

	as(as(0.5 * (lw + t(lw)), "symmetricMatrix"), "CsparseMatrix")
}
#}


