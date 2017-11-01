# Copyright 2001-10 by Roger Bivand, Markus Reder and Werner Mueller, 2015 Martin Gubri
#


nb2mat <- function(neighbours, glist=NULL, style="W", zero.policy=NULL)
{
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
	if(!inherits(neighbours, "nb")) stop("Not a neighbours list")
	listw <- nb2listw(neighbours, glist=glist, style=style,
		zero.policy=zero.policy)
	res <- listw2mat(listw)
	attr(res, "call") <- match.call()
	res
}

listw2mat <- function(listw) {
	n <- length(listw$neighbours)
	if (n < 1) stop("non-positive number of entities")
	cardnb <- card(listw$neighbours)
	if (any(is.na(unlist(listw$weights))))
		stop ("NAs in general weights list")
	res <- matrix(0, nrow=n, ncol=n)
	for (i in 1:n)
	    if (cardnb[i] > 0)
		res[i, listw$neighbours[[i]]] <- listw$weights[[i]]
	if (!is.null(attr(listw, "region.id")))
		row.names(res) <- attr(listw, "region.id")
	res
}

invIrM <- function(neighbours, rho, glist=NULL, style="W", method="solve", 
	feasible=NULL) {
	if(class(neighbours) != "nb") stop("Not a neighbours list")
	invIrW(nb2listw(neighbours, glist=glist, style=style), rho=rho, 
		method=method, feasible=feasible)
}

invIrW <- function(x, rho, method="solve", feasible=NULL) {
	if(inherits(x, "listw")) {
	  n <- length(x$neighbours)
	  V <- listw2mat(x)
	} else if (inherits(x, "Matrix") || inherits(x, "matrix")) {
	  if (method == "chol" && all(t(x) == x))
            stop("No Cholesky method for matrix or sparse matrix object")
          n <- dim(x)[1]
          V <- x
	} else stop("Not a weights list or a Sparse Matrix")
	if (is.null(feasible) || (is.logical(feasible) && !feasible)) {
		e <- eigen(V, only.values = TRUE)$values
		if (is.complex(e)) feasible <- 1/(range(Re(e)))
		else feasible <- 1/(range(e))
		if (rho <= feasible[1] || rho >= feasible[2])
			stop(paste("Rho", rho, "outside feasible range:",
                        paste(feasible, collapse=":")))
	}
	if (method == "chol"){
		if (x$style %in% c("W", "S") && !(can.be.simmed(x)))
			stop("Cholesky method requires symmetric weights")
		if (x$style %in% c("B", "C", "U") && 
			!(is.symmetric.glist(x$neighbours, x$weights)))
			stop("Cholesky method requires symmetric weights")
		if (x$style %in% c("W", "S")) {
			V <- listw2mat(listw2U(similar.listw(x)))
		}
		mat <- diag(n) - rho * V
		res <- chol2inv(chol(mat))
	} else if (method == "solve") {
		mat <- diag(n) - rho * V
		res <- solve(mat)
	} else stop("unknown method")
	attr(res, "call") <- match.call()
	res
}

powerWeights <- function(W, rho, order=250, X, tol=.Machine$double.eps^(3/5)) {
    timings <- list()
    .ptime_start <- proc.time()
    n <- dim(W)[1]
    dX <- dim(X)
    if (dX[1] == n) side <- "R"
    else if (dX[2] == n) side <- "L"
    else stop("W and X non-conformant")
    aW <- rho*W
    if (side == "R") last <- aW %*% X
    else last <- X %*% aW
    acc <- X + last
    conv <- FALSE
    iter <- 1
    series <- numeric(order)
    while (iter < order) {
        if (side == "R") {
            last <- aW %*% last
            acc <- acc + last
        } else {
            last <- last %*% aW
            acc <- acc + last
        }
# abs() added 2017-02-15, bug spotted by Yongwan Chun
        series[iter] <- mean(abs(as(last, "matrix")))
        if (series[iter] < tol) {
            conv <- TRUE
            break
        }
        iter <- iter+1
    }
    if (!conv) warning("not converged within order iterations")
    timings[["make_power_sum"]] <- proc.time() - .ptime_start
    attr(acc, "internal") <- list(series=series, order=order,
        tol=tol, iter=iter, conv=conv)
    attr(acc, "timings") <- do.call("rbind", timings)[, c(1, 3)]
    acc
}


mat2listw <- function(x, row.names=NULL, style="M") {
	if (!(is.matrix(x) || is(x, "sparseMatrix"))) stop("x is not a matrix")
	n <- nrow(x)
	if (n < 1) stop("non-positive number of entities")
	m <- ncol(x)
	if (n != m) stop("x must be a square matrix")
	if (any(x < 0)) stop("values in x cannot be negative")
	if (any(is.na(x))) stop("NA values in x not allowed")
    	if (!is.null(row.names)) {
		if(length(row.names) != n)
            		stop("row.names wrong length")
		if (length(unique(row.names)) != length(row.names))
	    		stop("non-unique row.names given")
    	}
    	if (is.null(row.names)) {
		if (!is.null(row.names(x))) {
			row.names <- row.names(x)
		} else {
			row.names <- as.character(1:n)
		}
	}
#	style <- "M"
        if (is(x, "sparseMatrix")) {
            xC <- as(x, "dgCMatrix")
            i <- slot(xC, "i")+1
            p <- slot(xC, "p")
            dp <- diff(p)
            rp <- rep(seq_along(dp), dp)
            df0 <- data.frame(from=i, to=rp, weights=slot(xC, "x"))
            o <- order(df0$from, df0$to)
            df <- df0[o,]
            class(df) <- c(class(df), "spatial.neighbour")
            attr(df, "region.id") <- row.names
            attr(df, "n") <- dim(xC)[1]
            res <- sn2listw(df)
            neighbours <- res$neighbours
            weights <- res$weights
        } else {
	    neighbours <- vector(mode="list", length=n)
	    weights <- vector(mode="list", length=n)
	    for (i in 1:n) {
		nbs  <- which(x[i,] > 0.0)
		if (length(nbs) > 0) {
			neighbours[[i]] <- nbs
			weights[[i]] <- as.double(x[i, nbs]) # Laurajean Lewis
		} else {
			neighbours[[i]] <- 0L
		}
	    }
        }
	attr(weights, "mode") <- "unknown" # Brian Rubineau
	class(neighbours) <- "nb"
	attr(neighbours, "region.id") <- row.names
 	attr(neighbours, "call") <- NA
        attr(neighbours, "sym") <- is.symmetric.nb(neighbours, 
		verbose=FALSE, force=TRUE)
	res <- list(style=style, neighbours=neighbours, weights=weights)
	class(res) <- c("listw", "nb")
	attr(res, "region.id") <- attr(neighbours, "region.id")
	attr(res, "call") <- match.call()
        if (style != "M") {
            res <- nb2listw(res$neighbours, glist=res$weights, style=style,
                zero.policy=TRUE)
        }
	res
}
