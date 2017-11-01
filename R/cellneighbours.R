# Copyright 2001 by Roger Bivand 
#

rookcell <- function(rowcol, nrow, ncol, torus=FALSE, rmin=1, cmin=1) {
	if (is.null(dim(rowcol))) rowcol <- t(as.matrix(rowcol))
        if(nrow(rowcol) != 1) stop("only single grid cell handled")
	row <- rowcol[1]
	col <- rowcol[2]
	if (torus) {
		y <- c(ifelse(col-1 < cmin, ncol, col-1), col, col,
			ifelse(col+1 > (ncol+(cmin-1)), cmin, col+1))
		x <- c(row, ifelse(row-1 < rmin, nrow, row-1),
			ifelse(row+1 > (nrow+(rmin-1)), rmin, row+1), row)
	} else {
		y <- c(ifelse(col-1 < cmin, NA, col-1), col, col,
			ifelse(col+1 > (ncol+(cmin-1)), NA, col+1))
		x <- c(row, ifelse(row-1 < rmin, NA, row-1),
			ifelse(row+1 > (nrow+(rmin-1)), NA, row+1), row)
	}
	res <- as.data.frame(list(row=x, col=y))
	res <- na.omit(res)
	res <- as.matrix(res)
	rownames(res) <- NULL
	attr(res, "coords") <- c(col, row)
	res
}

queencell <- function(rowcol, nrow, ncol, torus=FALSE, rmin=1, cmin=1) {
	if (is.null(dim(rowcol))) rowcol <- t(as.matrix(rowcol))
        if(nrow(rowcol) != 1) stop("only single grid cell handled")
 	row <- rowcol[1]
	col <- rowcol[2]
	if (torus) {
		y <- c(rep(ifelse(col-1 < cmin, ncol, col-1), 3), col, col,
			rep(ifelse(col+1 > (ncol+(cmin-1)), cmin, col+1), 3))
		x <- integer(8)
		x[c(1,4,6)] <- rep(ifelse(row+1 > (nrow+(rmin-1)),
			rmin, row+1), 3)
		x[c(2,7)] <- rep(row, 2)
		x[c(3,5,8)] <- rep(ifelse(row-1 < rmin, nrow, row-1), 3)
	} else {
		y <- c(rep(ifelse(col-1 < cmin, NA, col-1), 3), col, col,
			rep(ifelse(col+1 > (ncol+(cmin-1)), NA, col+1), 3))
		x <- integer(8)
		x[c(1,4,6)] <- rep(ifelse(row+1 > (nrow+(rmin-1)),
			NA, row+1), 3)
		x[c(2,7)] <- rep(row, 2)
		x[c(3,5,8)] <- rep(ifelse(row-1 < rmin, NA, row-1), 3)
	}
	res <- as.data.frame(list(row=x, col=y))
	res <- na.omit(res)
	res <- as.matrix(res)
	rownames(res) <- NULL
	attr(res, "coords") <- c(col, row)
	res
}

mrc2vi <- function(rowcol, nrow, ncol) {
	i <- ((rowcol[,2]-1) * nrow) + rowcol[,1]
	if (i > nrow*ncol || i < 1) stop("row or column out of range")
	as.integer(i)
}

vi2mrc <- function(i, nrow, ncol) {
	col <- ceiling(i/nrow)
	tmp <- i %% nrow
	row <- ifelse(tmp == 0, nrow, tmp)
	if (row < 1 || row > nrow) stop("i out of range")
	if (col < 1 || col > ncol) stop("i out of range")
	res <- cbind(row, col)
	res
}

cell2nb <- function(nrow, ncol, type="rook", torus=FALSE) {
	nrow <- as.integer(nrow)
	if (nrow < 1) stop("nrow nonpositive")
	ncol <- as.integer(ncol)
	if (ncol < 1) stop("nrow nonpositive")
	xcell <- NULL
	if (type == "rook") xcell <- rookcell
	if (type == "queen") xcell <- queencell
	if (is.null(xcell))
		stop(paste(type, ": no such cell function", sep=""))
	n <- nrow * ncol
	if (n < 0) stop("non-positive number of cells")
	res <- vector(mode="list", length=n)
	rownames <- character(n)
	for (i in 1:n) {
		res[[i]] <- sort(mrc2vi(xcell(vi2mrc(i, nrow, ncol),
			nrow, ncol, torus), nrow, ncol))
		rownames[i] <- paste(vi2mrc(i, nrow, ncol), collapse=":")
	}
	class(res) <- "nb"
	attr(res, "call") <- match.call()
	attr(res, "region.id") <- rownames
	attr(res, "cell") <- TRUE
	attr(res, type) <- TRUE
	if (torus) attr(res, "torus") <- TRUE
	res <- sym.attr.nb(res)
	res
}


