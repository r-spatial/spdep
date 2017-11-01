# Copyright 2003-2015 by Roger Bivand 

set.listw_is_CsparseMatrix_Option <- function(check) {
	if (!is.logical(check)) stop ("logical argument required")
	res <- get("listw_is_CsparseMatrix", envir = .spdepOptions)
	assign("listw_is_CsparseMatrix", check, envir = .spdepOptions)
	res
}

get.listw_is_CsparseMatrix_Option <- function() {
	get("listw_is_CsparseMatrix", envir = .spdepOptions)
}

set.spChkOption <- function(check) {
	if (!is.logical(check)) stop ("logical argument required")
	res <- get("spChkID", envir = .spdepOptions)
	assign("spChkID", check, envir = .spdepOptions)
	res
}

get.spChkOption <- function() {
	get("spChkID", envir = .spdepOptions)
}

set.VerboseOption <- function(check) {
	if (!is.logical(check)) stop ("logical argument required")
	res <- get("verbose", envir = .spdepOptions)
	assign("verbose", check, envir = .spdepOptions)
	res
}

get.VerboseOption <- function() {
	get("verbose", envir = .spdepOptions)
}

set.ZeroPolicyOption <- function(check) {
	if (!is.logical(check)) stop ("logical argument required")
	res <- get("zeroPolicy", envir = .spdepOptions)
	assign("zeroPolicy", check, envir = .spdepOptions)
	res
}

get.ZeroPolicyOption <- function() {
	get("zeroPolicy", envir = .spdepOptions)
}

set.ClusterOption <- function(cl) {
	if (!is.null(cl)) {
            if (!inherits(cl, "cluster")) stop ("cluster required")
        }
	assign("cluster", cl, envir = .spdepOptions)
        invisible(NULL)
}

get.ClusterOption  <- function() {
	get("cluster", envir = .spdepOptions)
}

set.mcOption <- function(value) {
        stopifnot(is.logical(value))
        stopifnot(length(value) == 1)
	res <- get("mc", envir = .spdepOptions)
        if (.Platform$OS.type == "windows") {
            if (value) warning("multicore not available on Windows")
        } else {
	    assign("mc", value, envir = .spdepOptions)
        }
	res
}

get.mcOption  <- function() {
	get("mc", envir = .spdepOptions)
}

set.coresOption <- function(value) {
	res <- get("cores", envir = .spdepOptions)
        if (is.null(value)) {
            assign("cores", value, envir = .spdepOptions)
        } else {
            stopifnot(is.integer(value))
            stopifnot(length(value) == 1)
            stopifnot(!is.na(value))
	    assign("cores", value, envir = .spdepOptions)
        }
	res
}

get.coresOption  <- function() {
	get("cores", envir = .spdepOptions)
}


chkIDs <- function (x, listw) 
{
    if (!is.array(x) & !is.data.frame(x)) {
        if (is.null(xn <- names(x))) 
            stop(paste(deparse(substitute(x)), "has no names"))
    }
    else {
        if (is.null(xn <- rownames(x))) 
            stop(paste(deparse(substitute(x)), "has no row names"))
    }
    if (!inherits(listw, "nb")) 
        stop(paste(deparse(substitute(listw)), "is not an listw  or nb object"))
    if (is.null(ln <- attr(listw, "region.id"))) 
        stop(paste(deparse(substitute(listw)), "has no region IDs"))
    if (length(ln) != length(xn)) 
        stop("objects of different length")
    res <- all(ln == xn)
    res
}

spNamedVec <- function(var, data) {
	if (!is.array(data) & !is.data.frame(data))
		stop(paste(deparse(substitute(data)),
			"not an array or data frame"))
	if (!is.character(var) & !is.numeric(var)) 
		stop("variable name wrong type") 
	res <- try(data[,var])
	if (inherits(res, "try-error")) 
		stop(paste(deparse(substitute(var)), "not found"))
	nms <- rownames(data)
	if (is.null(nms)) nms <- as.character(1:length(res))
	names(res) <- nms
	res
}
